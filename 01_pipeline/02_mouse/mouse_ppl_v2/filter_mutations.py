#!/usr/bin/env python3
import pandas as pd
import numpy as np
import gzip
import os
import sys
import argparse

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='从SNV文件中筛选germline和somatic mutation')
    
    parser.add_argument('-i', '--input', required=True, 
                       help='输入文件路径 (snv.tsv 或 snv.tsv.gz)')
    parser.add_argument('-o', '--output-dir', default=None,
                       help='输出目录路径 (默认为输入文件所在目录)')
    parser.add_argument('--germline-output', default='germline_mutations.tsv.gz',
                       help='germline mutations输出文件名 (默认: germline_mutations.tsv.gz)')
    parser.add_argument('--somatic-output', default='somatic_mutations.tsv.gz',
                       help='somatic mutations输出文件名 (默认: somatic_mutations.tsv.gz)')
    parser.add_argument('--report-output', default='mutation_filter_report.txt',
                       help='统计报告文件名 (默认: mutation_filter_report.txt)')
    parser.add_argument('--strand-min', type=float, default=0.3,
                       help='Strand_score最小值 (默认: 0.3)')
    parser.add_argument('--strand-max', type=float, default=0.7,
                       help='Strand_score最大值 (默认: 0.7)')
    parser.add_argument('--min-depth', type=int, default=10,
                       help='最小总reads数 (默认: 10)')
    parser.add_argument('--alt-ratio-threshold', type=float, default=0.9,
                       help='alt比例阈值，大于为germline，小于为somatic (默认: 0.9)')
    parser.add_argument('--no-compress', action='store_true',
                       help='不生成压缩文件')
    
    return parser.parse_args()

def filter_mutations(input_file, output_dir=None, germline_output='germline_mutations.tsv.gz',
                    somatic_output='somatic_mutations.tsv.gz', report_output='mutation_filter_report.txt',
                    strand_min=0.3, strand_max=0.7, min_depth=10, alt_ratio_threshold=0.9,
                    compress=True):
    """从snv文件中筛选germline和somatic mutation"""
    
    print(f"开始筛选突变")
    print(f"输入文件: {input_file}")
    print(f"输出目录: {output_dir}")
    print(f"筛选参数:")
    print(f"  Strand_score范围: {strand_min} - {strand_max}")
    print(f"  最小总reads数: {min_depth}")
    print(f"  alt比例阈值: {alt_ratio_threshold}")
    
    # 确定输出目录
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(input_file))
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取文件
    try:
        if input_file.endswith('.gz'):
            with gzip.open(input_file, 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"错误: 无法读取文件 {input_file}")
        print(f"错误信息: {e}")
        return None, None
    
    print(f"原始数据行数: {len(df)}")
    
    if df.empty:
        print("错误: 文件为空")
        return pd.DataFrame(), pd.DataFrame()
    
    # 显示列名确认
    print(f"列名: {list(df.columns)}")
    
    # 检查必需的列
    required_columns = ['Position', 'Ref', 'Alt', 'Ref_fw_total', 'Ref_rev_total',
                       'Alt_fw_total', 'Alt_rev_total', 'Strand_score', 'Mean_vaf']
    
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"错误: 缺少必需的列: {missing_columns}")
        return pd.DataFrame(), pd.DataFrame()
    
    # 计算必要的统计量
    df['total_reads'] = df['Ref_fw_total'] + df['Ref_rev_total'] + df['Alt_fw_total'] + df['Alt_rev_total']
    df['alt_reads'] = df['Alt_fw_total'] + df['Alt_rev_total']
    df['alt_ratio'] = df['alt_reads'] / df['total_reads']
    
    print(f"\n数据预览:")
    print(df[['Position', 'Ref', 'Alt', 'Strand_score', 'total_reads', 'alt_ratio', 'Mean_vaf']].head())
    
    # 筛选条件
    print(f"\n应用筛选条件:")
    
    # 条件1: Strand_score在指定范围内
    strand_condition = (df['Strand_score'] > strand_min) & (df['Strand_score'] < strand_max)
    print(f"条件1 (Strand_score在{strand_min}-{strand_max}之间): {strand_condition.sum()} 个突变")
    
    # 条件2: 总reads数大于阈值
    depth_condition = (df['total_reads'] > min_depth)
    print(f"条件2 (总reads > {min_depth}): {depth_condition.sum()} 个突变")
    
    # 条件3: alt比例 > 阈值 (germline)
    germline_condition = (df['alt_ratio'] > alt_ratio_threshold)
    print(f"条件3 (alt比例 > {alt_ratio_threshold}): {germline_condition.sum()} 个突变")
    
    # 条件4: alt比例 < 阈值 (somatic)
    somatic_condition = (df['alt_ratio'] < alt_ratio_threshold)
    print(f"条件4 (alt比例 < {alt_ratio_threshold}): {somatic_condition.sum()} 个突变")
    
    # 筛选germline mutation
    germline_df = df[strand_condition & depth_condition & germline_condition].copy()
    
    # 筛选somatic mutation
    somatic_df = df[strand_condition & depth_condition & somatic_condition].copy()
    
    print(f"\n筛选结果:")
    print(f"Germline mutations: {len(germline_df)}")
    print(f"Somatic mutations: {len(somatic_df)}")
    
    # 构建输出文件路径
    germline_path = os.path.join(output_dir, germline_output)
    somatic_path = os.path.join(output_dir, somatic_output)
    report_path = os.path.join(output_dir, report_output)
    
    # 保存结果
    if not germline_df.empty:
        if compress and germline_path.endswith('.gz'):
            germline_df.to_csv(germline_path, sep='\t', index=False, compression='gzip')
        else:
            # 如果指定了不压缩或者文件名没有.gz后缀，保存为普通文件
            output_path = germline_path.replace('.gz', '') if germline_path.endswith('.gz') else germline_path
            germline_df.to_csv(output_path, sep='\t', index=False)
            print(f"Germline mutations保存到: {output_path}")
    else:
        print("警告: 没有找到符合条件的germline mutations")
    
    if not somatic_df.empty:
        if compress and somatic_path.endswith('.gz'):
            somatic_df.to_csv(somatic_path, sep='\t', index=False, compression='gzip')
        else:
            # 如果指定了不压缩或者文件名没有.gz后缀，保存为普通文件
            output_path = somatic_path.replace('.gz', '') if somatic_path.endswith('.gz') else somatic_path
            somatic_df.to_csv(output_path, sep='\t', index=False)
            print(f"Somatic mutations保存到: {output_path}")
    else:
        print("警告: 没有找到符合条件的somatic mutations")
    
    # 生成详细统计报告
    with open(report_path, 'w') as f:
        f.write("=== 突变筛选报告 ===\n")
        f.write(f"生成时间: {pd.Timestamp.now()}\n")
        f.write(f"输入文件: {input_file}\n")
        f.write(f"输出目录: {output_dir}\n")
        f.write(f"\n筛选参数:\n")
        f.write(f"  Strand_score范围: {strand_min} - {strand_max}\n")
        f.write(f"  最小总reads数: {min_depth}\n")
        f.write(f"  alt比例阈值: {alt_ratio_threshold}\n")
        f.write(f"\n原始数据:\n")
        f.write(f"  总突变数: {len(df)}\n")
        f.write(f"  符合条件的突变数: {strand_condition.sum()}\n")
        f.write(f"  深度足够的突变数: {depth_condition.sum()}\n")
        f.write(f"\n筛选结果:\n")
        f.write(f"  Germline mutations: {len(germline_df)}\n")
        f.write(f"  Somatic mutations: {len(somatic_df)}\n")
        
        if not germline_df.empty:
            f.write(f"\nGermline统计:\n")
            f.write(f"  Strand_score范围: {germline_df['Strand_score'].min():.4f} - {germline_df['Strand_score'].max():.4f}\n")
            f.write(f"  平均Strand_score: {germline_df['Strand_score'].mean():.4f}\n")
            f.write(f"  alt比例范围: {germline_df['alt_ratio'].min():.4f} - {germline_df['alt_ratio'].max():.4f}\n")
            f.write(f"  平均alt比例: {germline_df['alt_ratio'].mean():.4f}\n")
            f.write(f"  平均总reads数: {germline_df['total_reads'].mean():.2f}\n")
            f.write(f"  平均VAF: {germline_df['Mean_vaf'].mean():.6f}\n")
            
            # 最常见的突变类型
            germline_df['mutation_type'] = germline_df['Ref'] + '>' + germline_df['Alt']
            mutation_counts = germline_df['mutation_type'].value_counts()
            f.write(f"\n  最常见的Germline突变类型 (前5):\n")
            for mutation, count in mutation_counts.head(5).items():
                f.write(f"    {mutation}: {count} 次\n")
        
        if not somatic_df.empty:
            f.write(f"\nSomatic统计:\n")
            f.write(f"  Strand_score范围: {somatic_df['Strand_score'].min():.4f} - {somatic_df['Strand_score'].max():.4f}\n")
            f.write(f"  平均Strand_score: {somatic_df['Strand_score'].mean():.4f}\n")
            f.write(f"  alt比例范围: {somatic_df['alt_ratio'].min():.4f} - {somatic_df['alt_ratio'].max():.4f}\n")
            f.write(f"  平均alt比例: {somatic_df['alt_ratio'].mean():.4f}\n")
            f.write(f"  平均总reads数: {somatic_df['total_reads'].mean():.2f}\n")
            f.write(f"  平均VAF: {somatic_df['Mean_vaf'].mean():.6f}\n")
            
            # 最常见的突变类型
            somatic_df['mutation_type'] = somatic_df['Ref'] + '>' + somatic_df['Alt']
            mutation_counts = somatic_df['mutation_type'].value_counts()
            f.write(f"\n  最常见的Somatic突变类型 (前5):\n")
            for mutation, count in mutation_counts.head(5).items():
                f.write(f"    {mutation}: {count} 次\n")
        
        # 保存前10个突变的详细信息
        if not germline_df.empty:
            f.write(f"\nGermline mutations (前10个):\n")
            f.write(germline_df.head(10).to_string(index=False))
            f.write("\n")
        
        if not somatic_df.empty:
            f.write(f"\nSomatic mutations (前10个):\n")
            f.write(somatic_df.head(10).to_string(index=False))
            f.write("\n")
    
    print(f"详细报告保存到: {report_path}")
    
    # 显示结果预览
    if not germline_df.empty:
        print(f"\nGermline mutations预览 (前5个):")
        print(germline_df[['Position', 'Ref', 'Alt', 'Strand_score', 'alt_ratio', 'Mean_vaf']].head().to_string(index=False))
    
    if not somatic_df.empty:
        print(f"\nSomatic mutations预览 (前5个):")
        print(somatic_df[['Position', 'Ref', 'Alt', 'Strand_score', 'alt_ratio', 'Mean_vaf']].head().to_string(index=False))
    
    return germline_df, somatic_df

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在 {args.input}")
        sys.exit(1)
    
    # 筛选突变
    germline_df, somatic_df = filter_mutations(
        input_file=args.input,
        output_dir=args.output_dir,
        germline_output=args.germline_output,
        somatic_output=args.somatic_output,
        report_output=args.report_output,
        strand_min=args.strand_min,
        strand_max=args.strand_max,
        min_depth=args.min_depth,
        alt_ratio_threshold=args.alt_ratio_threshold,
        compress=not args.no_compress
    )
    
    print(f"\n处理完成!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())