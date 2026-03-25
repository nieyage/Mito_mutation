#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import gzip
import argparse
from tqdm import tqdm
import sys
from collections import defaultdict
import multiprocessing as mp
from functools import partial
import tempfile
import pickle

def load_somatic_mutations(somatic_file):
    """加载somatic突变列表并转换为高效的查找结构"""
    print(f"加载somatic突变文件: {somatic_file}")
    
    try:
        if somatic_file.endswith('.gz'):
            with gzip.open(somatic_file, 'rt') as f:
                somatic_df = pd.read_csv(f, sep='\t')
        else:
            somatic_df = pd.read_csv(somatic_file, sep='\t')
    except Exception as e:
        print(f"错误: 无法读取somatic文件 {somatic_file}: {e}", file=sys.stderr)
        return None, None
    
    # 确保有必要的列
    required_cols = ['Position', 'Ref', 'Alt']
    missing_cols = [col for col in required_cols if col not in somatic_df.columns]
    
    if missing_cols:
        print(f"错误: somatic文件缺少列: {missing_cols}", file=sys.stderr)
        return None, None
    
    print(f"加载了 {len(somatic_df)} 个somatic突变")
    
    # 创建突变ID作为查找键
    somatic_df['mutation_id'] = somatic_df['Position'].astype(str) + '_' + somatic_df['Ref'] + '_' + somatic_df['Alt']
    
    # 创建两个查找字典：按位置和按突变ID
    mutations_by_pos = defaultdict(list)
    mutation_dict = {}
    
    for _, row in somatic_df.iterrows():
        pos = int(row['Position'])
        ref = row['Ref']
        alt = row['Alt']
        mutation_id = row['mutation_id']
        
        mutations_by_pos[pos].append((ref, alt, mutation_id))
        mutation_dict[mutation_id] = {'pos': pos, 'ref': ref, 'alt': alt}
    
    return mutations_by_pos, mutation_dict

def parse_single_snv_file(filepath, mutations_by_pos, mutation_dict):
    """解析单个snv文件并提取somatic突变信息"""
    try:
        barcode = os.path.basename(filepath).replace('.snv', '')
        
        # 读取文件
        df = pd.read_csv(filepath, sep='\t', comment='#')
        
        # 如果文件为空，返回所有突变为0的结果
        if df.empty:
            results = []
            for mutation_id, info in mutation_dict.items():
                results.append({
                    'barcode': barcode,
                    'pos': info['pos'],
                    'ref': info['ref'],
                    'alt': info['alt'],
                    'ref_count': 0,
                    'alt_count': 0,
                    'vaf': 0.0
                })
            return results
        
        # 创建结果的字典，初始化为0
        results_dict = {}
        for mutation_id, info in mutation_dict.items():
            results_dict[mutation_id] = {
                'barcode': barcode,
                'pos': info['pos'],
                'ref': info['ref'],
                'alt': info['alt'],
                'ref_count': 0,
                'alt_count': 0,
                'vaf': 0.0
            }
        
        # 只检查文件中存在的突变位置
        unique_positions = set(df['Position'].astype(int))
        positions_to_check = [pos for pos in unique_positions if pos in mutations_by_pos]
        
        # 如果没有需要检查的位置，直接返回所有为0的结果
        if not positions_to_check:
            return list(results_dict.values())
        
        # 批量处理需要检查的位置
        for pos in positions_to_check:
            pos_df = df[df['Position'] == pos]
            
            for ref, alt, mutation_id in mutations_by_pos[pos]:
                # 查找匹配的突变
                match = pos_df[(pos_df['Ref'] == ref) & (pos_df['VarAllele'] == alt)]
                
                if not match.empty:
                    row = match.iloc[0]
                    
                    # 获取reads计数
                    ref_count = int(row.get('Reads1', 0))
                    alt_count = int(row.get('Reads2', 0))
                    
                    # 处理VAF
                    vaf = 0.0
                    vaf_str = str(row.get('VarFreq', '0%')).replace('%', '')
                    try:
                        vaf = float(vaf_str) / 100.0
                    except:
                        vaf = 0.0
                    
                    # 更新结果
                    results_dict[mutation_id].update({
                        'ref_count': ref_count,
                        'alt_count': alt_count,
                        'vaf': vaf
                    })
        
        return list(results_dict.values())
        
    except Exception as e:
        print(f"解析文件 {filepath} 时出错: {e}", file=sys.stderr)
        # 出错时返回所有突变为0的结果
        results = []
        for mutation_id, info in mutation_dict.items():
            results.append({
                'barcode': os.path.basename(filepath).replace('.snv', ''),
                'pos': info['pos'],
                'ref': info['ref'],
                'alt': info['alt'],
                'ref_count': 0,
                'alt_count': 0,
                'vaf': 0.0
            })
        return results

def process_snv_files_parallel(snv_files, mutations_by_pos, mutation_dict, num_workers=None):
    """并行处理多个snv文件"""
    if num_workers is None:
        num_workers = max(1, mp.cpu_count() - 1)
    
    print(f"使用 {num_workers} 个进程并行处理 {len(snv_files)} 个文件...")
    
    # 创建部分函数，传入共享参数
    parse_func = partial(parse_single_snv_file, 
                        mutations_by_pos=mutations_by_pos, 
                        mutation_dict=mutation_dict)
    
    # 使用进程池并行处理
    all_results = []
    with mp.Pool(processes=num_workers) as pool:
        # 使用imap_unordered以获得更好的进度显示
        for results in tqdm(pool.imap_unordered(parse_func, snv_files), 
                          total=len(snv_files), 
                          desc="并行处理细胞"):
            all_results.extend(results)
    
    return all_results

def create_variant_sparse_matrix_batch(somatic_file, snv_dir, output_file, chromosome='chrM', 
                                      sample_cells=0, num_workers=None, batch_size=1000):
    """批量处理版本：创建variant稀疏矩阵"""
    
    print(f"开始创建variant稀疏矩阵...")
    print(f"Somatic文件: {somatic_file}")
    print(f"SNV目录: {snv_dir}")
    print(f"输出文件: {output_file}")
    print(f"染色体: {chromosome}")
    
    # 加载somatic突变
    mutations_by_pos, mutation_dict = load_somatic_mutations(somatic_file)
    if mutations_by_pos is None or mutation_dict is None:
        return False
    
    # 获取所有snv文件
    snv_files = glob.glob(os.path.join(snv_dir, "*.snv"))
    print(f"找到 {len(snv_files)} 个snv文件")
    
    if not snv_files:
        print("错误: 未找到snv文件", file=sys.stderr)
        return False
    
    # 如果需要抽样
    if sample_cells > 0 and sample_cells < len(snv_files):
        import random
        snv_files = random.sample(snv_files, sample_cells)
        print(f"测试模式: 只处理 {len(snv_files)} 个细胞")
    
    # 如果突变数量太多，考虑分批处理
    total_mutations = len(mutation_dict)
    print(f"需要跟踪 {total_mutations} 个somatic突变")
    
    # 创建输出文件
    compress = output_file.endswith('.gz')
    
    # 分批处理文件
    total_variants_written = 0
    num_batches = (len(snv_files) + batch_size - 1) // batch_size
    
    print(f"将 {len(snv_files)} 个文件分成 {num_batches} 批处理，每批最多 {batch_size} 个文件")
    
    with gzip.open(output_file, 'wt') if compress else open(output_file, 'w') as out_fh:
        # 写入表头
        out_fh.write("cell\tchrom\tpos\tref_base\talt_base\tref_count\talt_count\tvaf\n")
        
        # 分批处理文件
        for batch_idx in range(num_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(snv_files))
            batch_files = snv_files[start_idx:end_idx]
            
            print(f"\n处理批次 {batch_idx + 1}/{num_batches} ({len(batch_files)} 个文件)...")
            
            # 处理当前批次的文件
            batch_results = []
            if num_workers and num_workers > 1 and len(batch_files) > 1:
                batch_results = process_snv_files_parallel(batch_files, mutations_by_pos, 
                                                          mutation_dict, num_workers)
            else:
                # 单进程处理
                for snv_file in tqdm(batch_files, desc="处理文件"):
                    results = parse_single_snv_file(snv_file, mutations_by_pos, mutation_dict)
                    batch_results.extend(results)
            
            # 写入当前批次的结果
            for result in batch_results:
                out_fh.write(f"{result['barcode']}\t{chromosome}\t{result['pos']}\t"
                           f"{result['ref']}\t{result['alt']}\t"
                           f"{result['ref_count']}\t{result['alt_count']}\t{result['vaf']:.6f}\n")
                total_variants_written += 1
            
            # 清理内存
            del batch_results
    
    # 验证输出文件
    if os.path.exists(output_file):
        file_size = os.path.getsize(output_file) / (1024*1024*1024)  # GB
        if compress:
            # 压缩文件实际大小会小很多
            file_size = max(0.001, file_size * 0.1)  # 粗略估计
        
        print(f"\n处理完成!")
        print(f"总处理细胞数: {len(snv_files)}")
        print(f"总变异记录数: {total_variants_written}")
        print(f"预期记录数: {len(snv_files) * total_mutations}")
        print(f"输出文件大小: {file_size:.2f} GB")
        
        # 验证记录数
        expected_records = len(snv_files) * total_mutations
        if total_variants_written == expected_records:
            print("✓ 所有变异记录都已正确生成")
        else:
            print(f"⚠ 记录数不匹配: 生成 {total_variants_written}, 预期 {expected_records}")
        
        return True
    else:
        print("错误: 输出文件未生成", file=sys.stderr)
        return False

def create_variant_sparse_matrix_memory_efficient(somatic_file, snv_dir, output_file, 
                                                 chromosome='chrM', sample_cells=0):
    """内存高效版本：使用SQLite或分块处理"""
    
    print(f"使用内存高效模式创建variant稀疏矩阵...")
    
    # 加载somatic突变
    mutations_by_pos, mutation_dict = load_somatic_mutations(somatic_file)
    if mutations_by_pos is None or mutation_dict is None:
        return False
    
    # 获取所有snv文件
    snv_files = glob.glob(os.path.join(snv_dir, "*.snv"))
    print(f"找到 {len(snv_files)} 个snv文件")
    
    if not snv_files:
        print("错误: 未找到snv文件", file=sys.stderr)
        return False
    
    # 如果需要抽样
    if sample_cells > 0 and sample_cells < len(snv_files):
        import random
        snv_files = random.sample(snv_files, sample_cells)
        print(f"测试模式: 只处理 {len(snv_files)} 个细胞")
    
    total_mutations = len(mutation_dict)
    print(f"需要跟踪 {total_mutations} 个somatic突变")
    
    # 创建输出文件
    compress = output_file.endswith('.gz')
    
    # 预先创建所有突变的基础记录（用于填充缺失值）
    base_mutations = []
    mutation_ids = []
    for mutation_id, info in mutation_dict.items():
        base_mutations.append((info['pos'], info['ref'], info['alt']))
        mutation_ids.append(mutation_id)
    
    # 处理每个细胞文件
    total_variants_written = 0
    
    with (gzip.open(output_file, 'wt') if compress else open(output_file, 'w')) as out_fh:
        # 写入表头
        out_fh.write("cell\tchrom\tpos\tref_base\talt_base\tref_count\talt_count\tvaf\n")
        
        # 处理每个细胞
        for snv_file in tqdm(snv_files, desc="处理细胞"):
            barcode = os.path.basename(snv_file).replace('.snv', '')
            
            try:
                # 读取当前细胞的文件
                df = pd.read_csv(snv_file, sep='\t', comment='#')
                
                # 创建突变查找字典
                cell_mutations = {}
                if not df.empty:
                    # 处理文件中存在的突变
                    for _, row in df.iterrows():
                        pos = int(row['Position'])
                        ref = row['Ref']
                        alt = row['VarAllele']
                        
                        # 检查是否是somatic突变
                        mutation_id = f"{pos}_{ref}_{alt}"
                        if mutation_id in mutation_dict:
                            ref_count = int(row.get('Reads1', 0))
                            alt_count = int(row.get('Reads2', 0))
                            
                            # 处理VAF
                            vaf = 0.0
                            vaf_str = str(row.get('VarFreq', '0%')).replace('%', '')
                            try:
                                vaf = float(vaf_str) / 100.0
                            except:
                                vaf = 0.0
                            
                            cell_mutations[mutation_id] = (ref_count, alt_count, vaf)
                
                # 写入所有突变（存在的和缺失的）
                for i, mutation_id in enumerate(mutation_ids):
                    pos, ref, alt = base_mutations[i]
                    
                    if mutation_id in cell_mutations:
                        ref_count, alt_count, vaf = cell_mutations[mutation_id]
                    else:
                        ref_count, alt_count, vaf = 0, 0, 0.0
                    
                    out_fh.write(f"{barcode}\t{chromosome}\t{pos}\t{ref}\t{alt}\t"
                               f"{ref_count}\t{alt_count}\t{vaf:.6f}\n")
                    total_variants_written += 1
                    
            except Exception as e:
                print(f"处理文件 {snv_file} 时出错，使用默认值填充: {e}", file=sys.stderr)
                # 出错时填充所有突变为0
                for pos, ref, alt in base_mutations:
                    out_fh.write(f"{barcode}\t{chromosome}\t{pos}\t{ref}\t{alt}\t0\t0\t0.000000\n")
                    total_variants_written += 1
    
    print(f"\n处理完成!")
    print(f"总处理细胞数: {len(snv_files)}")
    print(f"总变异记录数: {total_variants_written}")
    print(f"预期记录数: {len(snv_files) * total_mutations}")
    
    return True

def main():
    parser = argparse.ArgumentParser(description='高效创建variant稀疏矩阵')
    parser.add_argument('-s', '--somatic-file', required=True,
                       help='somatic突变文件路径')
    parser.add_argument('-d', '--snv-dir', required=True,
                       help='包含所有细胞SNV文件的目录')
    parser.add_argument('-o', '--output', default='variant_sparse_matrix.tsv.gz',
                       help='输出文件路径 (默认: variant_sparse_matrix.tsv.gz)')
    parser.add_argument('--chromosome', default='chrM',
                       help='染色体名称 (默认: chrM)')
    parser.add_argument('--sample', type=int, default=0,
                       help='只处理指定数量的细胞进行测试 (0表示处理所有)')
    parser.add_argument('--workers', type=int, default=None,
                       help='并行处理的工作进程数 (默认: CPU核心数-1)')
    parser.add_argument('--batch-size', type=int, default=1000,
                       help='每批处理的文件数 (默认: 1000)')
    parser.add_argument('--memory-efficient', action='store_true',
                       help='使用内存高效模式 (处理大量突变时推荐)')
    parser.add_argument('--skip-preview', action='store_true',
                       help='跳过输出文件预览')
    
    args = parser.parse_args()
    
    # 创建输出目录
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 记录开始时间
    import time
    start_time = time.time()
    
    # 选择处理模式
    if args.memory_efficient:
        success = create_variant_sparse_matrix_memory_efficient(
            somatic_file=args.somatic_file,
            snv_dir=args.snv_dir,
            output_file=args.output,
            chromosome=args.chromosome,
            sample_cells=args.sample
        )
    else:
        success = create_variant_sparse_matrix_batch(
            somatic_file=args.somatic_file,
            snv_dir=args.snv_dir,
            output_file=args.output,
            chromosome=args.chromosome,
            sample_cells=args.sample,
            num_workers=args.workers,
            batch_size=args.batch_size
        )
    
    # 计算运行时间
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    if success:
        print(f"\n总运行时间: {elapsed_time:.2f} 秒")
        
        # 显示结果预览
        if not args.skip_preview and os.path.exists(args.output):
            print(f"\n输出文件预览 (前10行):")
            try:
                if args.output.endswith('.gz'):
                    with gzip.open(args.output, 'rt') as f:
                        lines = []
                        for i in range(10):
                            line = f.readline()
                            if line and i < 10:
                                lines.append(line.strip())
                        print('\n'.join(lines))
                else:
                    with open(args.output, 'r') as f:
                        lines = []
                        for i in range(10):
                            line = f.readline()
                            if line and i < 10:
                                lines.append(line.strip())
                        print('\n'.join(lines))
            except Exception as e:
                print(f"无法读取输出文件预览: {e}")
        
        # 生成统计报告
        stats_file = args.output.replace('.tsv', '_stats.txt').replace('.gz', '')
        try:
            with open(stats_file, 'w') as f:
                f.write("=== Variant Sparse Matrix 生成统计 ===\n")
                f.write(f"生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"输入目录: {args.snv_dir}\n")
                f.write(f"Somatic文件: {args.somatic_file}\n")
                f.write(f"输出文件: {args.output}\n")
                f.write(f"染色体: {args.chromosome}\n")
                f.write(f"运行时间: {elapsed_time:.2f} 秒\n")
                
                # 尝试读取输出文件获取更多统计信息
                if os.path.exists(args.output):
                    # 快速统计行数（不包括表头）
                    cmd = f"zcat {args.output} | wc -l" if args.output.endswith('.gz') else f"wc -l < {args.output}"
                    import subprocess
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if result.returncode == 0:
                        total_lines = int(result.stdout.strip())
                        f.write(f"总行数: {total_lines}\n")
        
            print(f"统计报告保存至: {stats_file}")
        except Exception as e:
            print(f"无法生成统计报告: {e}")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
