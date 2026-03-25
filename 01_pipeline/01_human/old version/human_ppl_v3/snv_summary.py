#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import gzip

def parse_snv_file(filepath):
    """解析单个snv文件"""
    barcode = os.path.basename(filepath).replace('.snv', '')
    
    try:
        # 读取varscan输出文件
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

def main():
    input_dir = "/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell"
    output_file = f"./snv.tsv"
    min_depth = 2
    
    print("开始处理所有细胞的SNV结果...")
    
    # 获取所有snv文件
    snv_files = glob.glob(f"{input_dir}/*.snv")
    print(f"找到 {len(snv_files)} 个snv文件")
    
    if not snv_files:
        print("错误: 未找到snv文件")
        return 1
    
    # 步骤1: 解析所有文件
    print("步骤1: 解析所有snv文件...")
    all_variants = []
    
    for i, snv_file in enumerate(snv_files):
        df = parse_snv_file(snv_file)
        if not df.empty:
            all_variants.append(df)
        
        # 显示进度
        if (i + 1) % 1000 == 0:
            print(f"已解析 {i + 1}/{len(snv_files)} 个文件")
    
    # 合并所有变异数据
    if not all_variants:
        print("错误: 未解析到任何变异数据")
        return 1
    
    combined_df = pd.concat(all_variants, ignore_index=True)
    print(f"总变异记录数: {len(combined_df)}")
    
    # 步骤2: 确定每个位置的alt_base
    print("步骤2: 确定每个位置的alt_base...")
    
    # 按位置和ref分组，找到最丰富的alt碱基
    mutation_stats = []
    
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
        
        # 找出最丰富的alt碱基
        alt_counts = pos_variants.groupby('alt').apply(
            lambda x: (x['alt_fw'] + x['alt_rev']).sum()
        )
        
        if not alt_counts.empty:
            dominant_alt = alt_counts.idxmax()
            total_count = alt_counts.max()
            
            mutation_stats.append({
                'position': position,
                'ref': ref_base,
                'alt': dominant_alt,
                'total_count': total_count
            })
    
    print(f"确定的突变数量: {len(mutation_stats)}")
    
    # 步骤3: 对每个突变进行计算
    print("步骤3: 计算突变统计量...")
    
    results = []
    
    for mutation in mutation_stats:
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
        
        n_cells = len(mutation_cells)
        conf_cell = 0
        vaf_list = []
        vaf_pos_cells = 0
        
        # 初始化总计
        ref_fw_total = 0
        ref_rev_total = 0
        alt_fw_total = 0
        alt_rev_total = 0
        
        # 处理每个细胞
        for _, row in mutation_cells.iterrows():
            ref_fw = int(row['ref_fw'])
            ref_rev = int(row['ref_rev'])
            alt_fw = int(row['alt_fw'])
            alt_rev = int(row['alt_rev'])
            vaf = float(row['vaf'])
            
            # 累加总reads数
            ref_fw_total += ref_fw
            ref_rev_total += ref_rev
            alt_fw_total += alt_fw
            alt_rev_total += alt_rev
            
            # 检查是否为confident cell
            if alt_fw >= min_depth and alt_rev >= min_depth:
                conf_cell += 1
            
            # 收集VAF
            vaf_list.append(vaf)
            
            # 检查是否有变异
            if vaf > 0:
                vaf_pos_cells += 1
        
        # 计算统计量
        if vaf_list:
            vaf_array = np.array(vaf_list)
            mean_vaf = vaf_array.mean()
            var_vaf = vaf_array.var()
            lis = mean_vaf * (1 + var_vaf)
        else:
            mean_vaf = 0
            var_vaf = 0
            lis = 0
        
        # 计算百分比
        pct_conf = conf_cell / n_cells if n_cells > 0 else 0
        pct_vaf_pos = vaf_pos_cells / n_cells if n_cells > 0 else 0
        
        # 计算链偏好性分数
        strand_score = 0
        if alt_fw_total + alt_rev_total > 0:
            strand_score = alt_fw_total / (alt_fw_total + alt_rev_total)
        
        if conf_cell > 0:
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
                'Pct_conf': round(pct_conf, 4),
                'Pct_vaf_pos': round(pct_vaf_pos, 4),
                'N_cells': n_cells,
                'Conf_cells': conf_cell,
                'VAF_pos_cells': vaf_pos_cells
            })
    
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('Position')

    # 保存结果
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"结果已保存到: {output_file}")
    
    # 显示统计信息
    print(f"\n=== 统计报告 ===")
    print(f"总突变数: {len(results_df)}")
    print(f"平均VAF: {results_df['Mean_vaf'].mean():.6f}")
    print(f"平均confident细胞比例: {results_df['Pct_conf'].mean():.4f}")
    print(f"平均VAF阳性细胞比例: {results_df['Pct_vaf_pos'].mean():.4f}")
    # 可选：保存为压缩格式
    results_df.to_csv(f"{output_file}.gz", sep='\t', index=False, compression='gzip')
    print(f"压缩版本已保存到: {output_file}.gz")
    
    return 0

if __name__ == "__main__":
    exit(main())