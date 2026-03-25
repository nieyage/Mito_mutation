#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import gzip

def parse_snv_file(filepath):
    barcode = os.path.basename(filepath).replace('.snv', '')
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#')
        column_map = {
            'Chrom': 'chrom',
            'Position': 'position',
            'Ref': 'ref',
            'VarAllele': 'alt',
            'Reads1': 'ref_total',
            'Reads2': 'alt_total',
            'VarFreq': 'vaf',
            'Reads1Plus': 'ref_fw',
            'Reads1Minus': 'ref_rev',
            'Reads2Plus': 'alt_fw',
            'Reads2Minus': 'alt_rev'
        }
        
        df = df.rename(columns=column_map)
        df['barcode'] = barcode
        if 'vaf' in df.columns:
            df['vaf'] = df['vaf'].astype(str).str.replace('%', '', regex=False)
            df['vaf'] = pd.to_numeric(df['vaf'], errors='coerce') / 100.0
        required_cols = ['barcode', 'position', 'ref', 'alt', 
                        'ref_fw', 'ref_rev', 'alt_fw', 'alt_rev', 'vaf']

        available_cols = [col for col in required_cols if col in df.columns]
        df = df[available_cols]
        
        return df
        
    except Exception as e:
        print(f"解析文件 {filepath} 时出错: {e}")
        return pd.DataFrame()

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='处理所有细胞的SNV结果，计算突变统计量',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python %(prog)s -i /path/to/snv/files -o ./snv_results.tsv
  python %(prog)s -i /path/to/snv/files -o ./snv_results.tsv --min-depth 3
  python %(prog)s -i /path/to/snv/files -o ./snv_results.tsv --no-compress
        '''
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        type=str,
        required=True,
        help='输入目录，包含所有.snv文件'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='./snv.tsv',
        help='输出文件路径 (默认: ./snv.tsv)'
    )
    
    parser.add_argument(
        '--min-depth',
        type=int,
        default=2,
        help='confident cell的最小深度要求 (默认: 2)'
    )
    
    parser.add_argument(
        '--pattern',
        type=str,
        default='*.snv',
        help='文件匹配模式 (默认: *.snv)'
    )
    
    parser.add_argument(
        '--no-compress',
        action='store_true',
        help='不生成压缩文件'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='显示详细进度信息'
    )
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_dir = args.input_dir
    output_file = args.output
    min_depth = args.min_depth
    pattern = args.pattern
    compress = not args.no_compress
    verbose = args.verbose
    
    print("开始处理所有细胞的SNV结果...")
    print(f"输入目录: {input_dir}")
    print(f"输出文件: {output_file}")
    print(f"最小深度: {min_depth}")
    print(f"文件模式: {pattern}")
    print(f"生成压缩文件: {'是' if compress else '否'}")
    
    if not os.path.exists(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        return 1
    
    snv_files = glob.glob(os.path.join(input_dir, pattern))
    print(f"找到 {len(snv_files)} 个snv文件")
    
    if not snv_files:
        print(f"错误: 在 {input_dir} 中未找到匹配 {pattern} 的文件")
        return 1
    
    print("步骤1: 解析所有snv文件...")
    all_variants = []
    all_cell_barcodes = set()

    for i, snv_file in enumerate(snv_files):
        df = parse_snv_file(snv_file)
        if not df.empty:
            all_variants.append(df)
            all_cell_barcodes.add(os.path.basename(snv_file).replace('.snv', ''))
        
        # 显示进度
        if verbose and (i + 1) % 1000 == 0:
            print(f"已解析 {i + 1}/{len(snv_files)} 个文件")
    
    # 合并所有变异数据
    if not all_variants:
        print("错误: 未解析到任何变异数据")
        return 1
    
    combined_df = pd.concat(all_variants, ignore_index=True)
    print(f"总变异记录数: {len(combined_df)}")
    
    # 获取总细胞数
    total_cells = len(all_cell_barcodes)
    print(f"总细胞数: {total_cells}")
    
    # 步骤2: 确定所有非参考碱基作为候选突变
    print("步骤2: 收集所有非参考碱基作为候选突变...")
    
    # 按位置和ref分组，收集所有非参考alt碱基
    candidate_mutations = []
    
    # 获取所有唯一的位置和ref组合
    unique_positions = combined_df[['position', 'ref']].drop_duplicates()
    
    for idx, (position, ref_base) in unique_positions.iterrows():
        # 获取该位置的所有变异
        pos_variants = combined_df[
            (combined_df['position'] == position) & 
            (combined_df['ref'] == ref_base)
        ]
        
        if pos_variants.empty:
            continue
        
        # 收集所有非参考alt碱基
        unique_alts = pos_variants['alt'].unique()
        
        for alt_base in unique_alts:
            # 确保alt不是参考碱基
            if alt_base != ref_base:
                candidate_mutations.append({
                    'position': position,
                    'ref': ref_base,
                    'alt': alt_base
                })
    
    print(f"候选突变数量: {len(candidate_mutations)}")
    
    # 步骤3: 筛选突变 - 要求confident cell >= 1
    print("步骤3: 筛选突变 (confident cell >= 1)...")
    
    filtered_mutations = []
    
    for i, mutation in enumerate(candidate_mutations):
        position = mutation['position']
        ref_base = mutation['ref']
        alt_base = mutation['alt']
        
        # 获取该突变在所有细胞中的数据
        mutation_cells = combined_df[
            (combined_df['position'] == position) & 
            (combined_df['ref'] == ref_base) & 
            (combined_df['alt'] == alt_base)
        ]
        
        if mutation_cells.empty:
            continue
        
        # 计算confident cell数量 (alt_fw>=2 & alt_rev>=2)
        confident_cells = 0
        for _, row in mutation_cells.iterrows():
            alt_fw = int(row['alt_fw'])
            alt_rev = int(row['alt_rev'])
            if alt_fw >= min_depth and alt_rev >= min_depth:
                confident_cells += 1
        
        # 筛选条件：confident cell >= 1
        if confident_cells >= 1:
            filtered_mutations.append(mutation)
    
    print(f"筛选后突变数量: {len(filtered_mutations)}")
    
    if not filtered_mutations:
        print("错误: 没有突变满足筛选条件 (confident cell >= 1)")
        return 1
    
    # 步骤4: 对每个筛选后的突变计算统计量
    print("步骤4: 计算突变统计量...")
    
    results = []
    
    for i, mutation in enumerate(filtered_mutations):
        position = mutation['position']
        ref_base = mutation['ref']
        alt_base = mutation['alt']
        
        # 获取该突变在所有细胞中的数据
        mutation_cells = combined_df[
            (combined_df['position'] == position) & 
            (combined_df['ref'] == ref_base) & 
            (combined_df['alt'] == alt_base)
        ]
        
        if mutation_cells.empty:
            continue
        
        # 初始化总计（基于所有细胞）
        ref_fw_total = 0
        ref_rev_total = 0
        alt_fw_total = 0
        alt_rev_total = 0
        
        # 用于VAF统计（基于VAF>0的细胞）
        vaf_pos_cells_list = []
        vaf_pos_cells_count = 0
        confident_cells_in_vaf_pos = 0
        
        # 处理每个细胞
        for _, row in mutation_cells.iterrows():
            ref_fw = int(row['ref_fw'])
            ref_rev = int(row['ref_rev'])
            alt_fw = int(row['alt_fw'])
            alt_rev = int(row['alt_rev'])
            vaf = float(row['vaf'])
            
            # 累加总reads数（所有细胞）
            ref_fw_total += ref_fw
            ref_rev_total += ref_rev
            alt_fw_total += alt_fw
            alt_rev_total += alt_rev
            
            # VAF > 0 的细胞
            if vaf > 0:
                vaf_pos_cells_count += 1
                vaf_pos_cells_list.append(vaf)
                
                # 检查是否为confident cell (在VAF>0的细胞中)
                if alt_fw >= min_depth and alt_rev >= min_depth:
                    confident_cells_in_vaf_pos += 1
        
        # 计算统计量
        
        # 计算链偏好性分数（基于所有细胞的reads）
        strand_score = 0
        if alt_fw_total + alt_rev_total > 0:
            strand_score = alt_fw_total / (alt_fw_total + alt_rev_total)
        
        # 计算VAF统计量（基于VAF>0的细胞）
        if vaf_pos_cells_list:
            vaf_array = np.array(vaf_pos_cells_list)
            mean_vaf = vaf_array.mean()
            var_vaf = vaf_array.var()
            lis = mean_vaf / (1 + var_vaf) if var_vaf >= 0 else 0
        else:
            mean_vaf = 0
            var_vaf = 0
            lis = 0
        
        # 计算比例
        # Pct_conf: VAF>0的细胞中confident cell的比例
        pct_conf = confident_cells_in_vaf_pos / vaf_pos_cells_count if vaf_pos_cells_count > 0 else 0
        
        # Pct_vaf_pos: 所有细胞中VAF>0的细胞比例
        pct_vaf_pos = vaf_pos_cells_count / total_cells if total_cells > 0 else 0
        
        results.append({
            'Position': position,
            'Ref': ref_base,
            'Alt': alt_base,
            'Ref_fw_total': ref_fw_total,
            'Ref_rev_total': ref_rev_total,
            'Alt_fw_total': alt_fw_total,
            'Alt_rev_total': alt_rev_total,
            'Strand_score': round(strand_score, 4),
            'Mean_vaf': round(mean_vaf, 6),
            'Var_vaf': round(var_vaf, 6),
            'Lis': round(lis, 6),
            'Pct_conf': round(pct_conf, 4),  # VAF>0的细胞中confident cell的比例
            'Pct_vaf_pos': round(pct_vaf_pos, 4),  # 所有细胞中VAF>0的细胞比例
            'Conf_cells_in_vaf_pos': confident_cells_in_vaf_pos,  # VAF>0细胞中的confident cell数
            'N_cells_vaf_pos': vaf_pos_cells_count,  # VAF>0的细胞数
            'Total_cells': total_cells  # 总细胞数
        })
        
        # 显示进度
        if verbose and (i + 1) % 1000 == 0:
            print(f"已处理 {i + 1}/{len(filtered_mutations)} 个突变")
    
    if not results:
        print("错误: 未计算出任何突变结果")
        return 1
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Position')
    
    # 保存结果
    print(f"保存结果到: {output_file}")
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # 可选：保存为压缩格式
    if compress:
        compressed_file = f"{output_file}.gz"
        print(f"保存压缩版本到: {compressed_file}")
        results_df.to_csv(compressed_file, sep='\t', index=False, compression='gzip')
    
    # 显示统计信息
    print(f"\n=== 统计报告 ===")
    print(f"总突变数: {len(results_df)}")
    print(f"平均VAF: {results_df['Mean_vaf'].mean():.6f}")
    print(f"平均Pct_conf (VAF>0细胞中confident细胞比例): {results_df['Pct_conf'].mean():.4f}")
    print(f"平均Pct_vaf_pos (所有细胞中VAF>0比例): {results_df['Pct_vaf_pos'].mean():.4f}")
    print(f"平均VAF阳性细胞数: {results_df['N_cells_vaf_pos'].mean():.2f}")
    
    return 0

if __name__ == "__main__":
    exit(main())