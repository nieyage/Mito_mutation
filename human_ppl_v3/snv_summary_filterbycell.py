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
        
        # 首先应用列映射
        column_map = {
            'Chrom': 'chrom',
            'Position': 'position',
            'Ref': 'ref',
            'VarAllele': 'alt',
            'Reads1': 'Reads1',
            'Reads2': 'Reads2',
            'VarFreq': 'vaf',
            'Reads1Plus': 'Reads1Plus',
            'Reads1Minus': 'Reads1Minus',
            'Reads2Plus': 'Reads2Plus',
            'Reads2Minus': 'Reads2Minus'
        }
        
        df = df.rename(columns=column_map)
        
        # 筛选条件：只保留符合条件的行
        # 条件1: Reads2Plus > 1 & Reads2Minus > 1
        condition1 = (df['Reads2Plus'] > 1) & (df['Reads2Minus'] > 1)
        
        # 条件2: Reads2Plus/(Reads2Minus + Reads2Plus) < 0.7
        # 条件3: Reads2Plus/(Reads2Minus + Reads2Plus) > 0.3
        reads2_total = df['Reads2Plus'] + df['Reads2Minus']
        reads2_ratio = df['Reads2Plus'] / reads2_total.replace(0, np.nan)  # 避免除以0
        condition2 = reads2_ratio < 0.7
        condition3 = reads2_ratio > 0.3
        
        # 条件4: (Reads1 + Reads2) >= 10
        condition4 = (df['Reads1'] + df['Reads2']) >= 10
        
        # 条件5: Reads2/(Reads1 + Reads2) >= 0.10
        total_reads = df['Reads1'] + df['Reads2']
        var_ratio = df['Reads2'] / total_reads.replace(0, np.nan)  # 避免除以0
        condition5 = var_ratio >= 0.10
        
        # 应用所有条件
        filtered_df = df[condition1 & condition2 & condition3 & condition4 & condition5].copy()
        
        if filtered_df.empty:
            # print(f"文件 {barcode} 没有符合条件的变异")
            return pd.DataFrame()
        
        # 打印筛选统计（每100个文件显示一次）
        global file_counter
        if 'file_counter' not in globals():
            file_counter = 0
        
        file_counter += 1
        if file_counter % 1000 == 0:
            print(f"已处理 {file_counter} 个文件，当前文件 {barcode}: 原始行数={len(df)}, 筛选后行数={len(filtered_df)}")
        
        # 现在重命名列用于后续处理
        final_column_map = {
            'chrom': 'chrom',
            'position': 'position',
            'ref': 'ref',
            'alt': 'alt',
            'Reads1': 'ref_total',
            'Reads2': 'alt_total',
            'vaf': 'vaf',
            'Reads1Plus': 'ref_fw',
            'Reads1Minus': 'ref_rev',
            'Reads2Plus': 'alt_fw',
            'Reads2Minus': 'alt_rev'
        }
        
        filtered_df = filtered_df.rename(columns=final_column_map)
        filtered_df['barcode'] = barcode
        
        # 处理VAF百分比格式
        if 'vaf' in filtered_df.columns:
            filtered_df['vaf'] = filtered_df['vaf'].astype(str).str.replace('%', '', regex=False)
            filtered_df['vaf'] = pd.to_numeric(filtered_df['vaf'], errors='coerce') / 100.0
        
        # 确保数值列的数据类型正确
        numeric_cols = ['position', 'ref_fw', 'ref_rev', 'alt_fw', 'alt_rev']
        for col in numeric_cols:
            if col in filtered_df.columns:
                filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce').fillna(0).astype(int)
        
        # 确保vaf是浮点数
        if 'vaf' in filtered_df.columns:
            filtered_df['vaf'] = pd.to_numeric(filtered_df['vaf'], errors='coerce').fillna(0.0)
        
        # 选择需要的列
        required_cols = ['barcode', 'position', 'ref', 'alt', 
                        'ref_fw', 'ref_rev', 'alt_fw', 'alt_rev', 'vaf']
        
        available_cols = [col for col in required_cols if col in filtered_df.columns]
        filtered_df = filtered_df[available_cols]
        
        return filtered_df
        
    except Exception as e:
        print(f"解析文件 {filepath} 时出错: {e}")
        return pd.DataFrame()

def main():
    input_dir = "/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell"
    output_dir ="."
    output_file = f"{output_dir}/snv_filtered.tsv"
    min_depth = 2
    
    print("开始处理所有细胞的SNV结果...")
    print("应用筛选条件:")
    print("  1. Reads2Plus > 1 & Reads2Minus > 1")
    print("  2. Reads2Plus/(Reads2Minus+Reads2Plus) < 0.7")
    print("  3. Reads2Plus/(Reads2Minus+Reads2Plus) > 0.3")
    print("  4. (Reads1 + Reads2) >= 10")
    print("  5. Reads2/(Reads1 + Reads2) >= 0.10")
    
    # 获取所有snv文件
    snv_files = glob.glob(f"{input_dir}/*.snv")
    print(f"找到 {len(snv_files)} 个snv文件")
    
    if not snv_files:
        print("错误: 未找到snv文件")
        return 1
    
    # 步骤1: 解析所有文件
    print("步骤1: 解析并筛选所有snv文件...")
    all_variants = []
    total_original_rows = 0
    total_filtered_rows = 0
    
    for i, snv_file in enumerate(snv_files):
        df = parse_snv_file(snv_file)
        if not df.empty:
            all_variants.append(df)
            total_filtered_rows += len(df)
        
        # 显示进度
        if (i + 1) % 1000 == 0:
            print(f"已解析 {i + 1}/{len(snv_files)} 个文件")
    
    # 合并所有变异数据
    if not all_variants:
        print("错误: 未解析到任何符合条件的变异数据")
        return 1
    
    combined_df = pd.concat(all_variants, ignore_index=True)
    print(f"\n筛选统计:")
    print(f"总变异记录数 (筛选后): {len(combined_df)}")
    
    # 步骤2: 确定每个位置的alt_base
    print("\n步骤2: 确定每个位置的alt_base...")
    
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
    print("\n步骤3: 计算突变统计量...")
    
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
                'Conf_cells': conf_cell,
                'VAF_pos_cells': vaf_pos_cells
            })
    
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('Position')
        
        # 保存结果
        results_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n结果已保存到: {output_file}")
        
        # 显示统计信息
        print(f"\n=== 统计报告 ===")
        print(f"总突变数: {len(results_df)}")
        print(f"平均VAF: {results_df['Mean_vaf'].mean():.6f}")
        print(f"平均confident细胞比例: {results_df['Pct_conf'].mean():.4f}")
        print(f"平均VAF阳性细胞比例: {results_df['Pct_vaf_pos'].mean():.4f}")
        
        # 显示前10个突变
        print(f"\n前10个突变:")
        print(results_df.head(10).to_string(index=False))
        
        # 可选：保存为压缩格式
        results_df.to_csv(f"{output_file}.gz", sep='\t', index=False, compression='gzip')
        print(f"压缩版本已保存到: {output_file}.gz")
    else:
        print("警告: 没有找到Conf_cells > 0的突变")
        # 创建空的输出文件
        pd.DataFrame().to_csv(output_file, sep='\t', index=False)
    
    return 0

if __name__ == "__main__":
    exit(main())