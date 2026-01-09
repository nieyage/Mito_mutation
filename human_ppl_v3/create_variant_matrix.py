#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import gzip
import argparse
from tqdm import tqdm
import sys

def load_somatic_mutations(somatic_file):
    """加载somatic突变列表"""
    print(f"加载somatic突变文件: {somatic_file}")
    
    if somatic_file.endswith('.gz'):
        with gzip.open(somatic_file, 'rt') as f:
            somatic_df = pd.read_csv(f, sep='\t')
    else:
        somatic_df = pd.read_csv(somatic_file, sep='\t')
    
    # 确保有必要的列
    required_cols = ['Position', 'Ref', 'Alt']
    missing_cols = [col for col in required_cols if col not in somatic_df.columns]
    
    if missing_cols:
        print(f"错误: somatic文件缺少列: {missing_cols}", file=sys.stderr)
        return []
    
    # 转换为列表格式，提高查询效率
    somatic_mutations = []
    for _, row in somatic_df.iterrows():
        somatic_mutations.append({
            'position': row['Position'],
            'ref': row['Ref'],
            'alt': row['Alt']
        })
    
    print(f"加载了 {len(somatic_mutations)} 个somatic突变")
    return somatic_mutations

def parse_snv_file_fast(filepath, somatic_mutations):
    """快速解析单个snv文件并提取相关信息"""
    barcode = os.path.basename(filepath).replace('.snv', '')
    
    try:
        # 读取文件
        df = pd.read_csv(filepath, sep='\t', comment='#')
        
        # 为每个somatic突变查找信息
        results = []
        for mutation in somatic_mutations:
            pos = mutation['position']
            ref = mutation['ref']
            alt = mutation['alt']
            
            # 查找匹配的突变
            match = df[(df['Position'] == pos) & (df['Ref'] == ref) & (df['VarAllele'] == alt)]
            
            if not match.empty:
                row = match.iloc[0]
                ref_count = int(row['Reads1']) if 'Reads1' in df.columns else 0
                alt_count = int(row['Reads2']) if 'Reads2' in df.columns else 0
                
                # 处理VAF
                vaf = 0.0
                if 'VarFreq' in df.columns:
                    vaf_str = str(row['VarFreq'])
                    vaf_str = vaf_str.replace('%', '')
                    try:
                        vaf = float(vaf_str) / 100.0
                    except:
                        vaf = 0.0
            else:
                ref_count = 0
                alt_count = 0
                vaf = 0.0
            
            results.append((barcode, pos, ref, alt, ref_count, alt_count, vaf))
        
        return results
        
    except Exception as e:
        print(f"解析文件 {filepath} 时出错: {e}", file=sys.stderr)
        return []

def create_variant_sparse_matrix_optimized(somatic_file, snv_dir, output_file, chromosome='chrM', sample_cells=0):
    """优化版本：创建variant稀疏矩阵"""
    
    print(f"开始创建variant稀疏矩阵...")
    print(f"Somatic文件: {somatic_file}")
    print(f"SNV目录: {snv_dir}")
    print(f"输出文件: {output_file}")
    print(f"染色体: {chromosome}")
    
    # 加载somatic突变
    somatic_mutations = load_somatic_mutations(somatic_file)
    if not somatic_mutations:
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
    
    # 创建输出文件
    compress = output_file.endswith('.gz')
    if compress:
        out_fh = gzip.open(output_file, 'wt')
    else:
        out_fh = open(output_file, 'w')
    
    # 写入表头
    out_fh.write("cell\tchrom\tpos\tref_base\talt_base\tref_count\talt_count\tvaf\n")
    
    # 处理每个snv文件
    total_variants_written = 0
    
    for snv_file in tqdm(snv_files, desc="处理细胞"):
        # 解析文件并获取变异信息
        variants = parse_snv_file_fast(snv_file, somatic_mutations)
        
        for barcode, pos, ref, alt, ref_count, alt_count, vaf in variants:
            out_fh.write(f"{barcode}\t{chromosome}\t{pos}\t{ref}\t{alt}\t{ref_count}\t{alt_count}\t{vaf:.6f}\n")
            total_variants_written += 1
    
    out_fh.close()
    
    print(f"\n处理完成!")
    print(f"总处理细胞数: {len(snv_files)}")
    print(f"总变异记录数: {total_variants_written}")
    print(f"输出文件大小: {os.path.getsize(output_file) / (1024*1024):.2f} MB")
    
    return True

def main():
    parser = argparse.ArgumentParser(description='创建variant稀疏矩阵')
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
    
    args = parser.parse_args()
    
    # 创建输出目录
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 创建variant稀疏矩阵
    success = create_variant_sparse_matrix_optimized(
        somatic_file=args.somatic_file,
        snv_dir=args.snv_dir,
        output_file=args.output,
        chromosome=args.chromosome,
        sample_cells=args.sample
    )
    
    # 显示结果预览
    if success and os.path.exists(args.output):
        print(f"\n输出文件预览 (前10行):")
        try:
            if args.output.endswith('.gz'):
                with gzip.open(args.output, 'rt') as f:
                    lines = []
                    for i in range(10):
                        line = f.readline()
                        if line:
                            lines.append(line.strip())
                    print('\n'.join(lines))
            else:
                with open(args.output, 'r') as f:
                    lines = []
                    for i in range(10):
                        line = f.readline()
                        if line:
                            lines.append(line.strip())
                    print('\n'.join(lines))
        except Exception as e:
            print(f"无法读取输出文件预览: {e}")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())