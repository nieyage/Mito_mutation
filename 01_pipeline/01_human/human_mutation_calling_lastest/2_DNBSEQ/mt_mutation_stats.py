import pandas as pd
import numpy as np
import anndata as ad
import argparse
import os
from scipy.sparse import issparse

def get_layer_data(adata, layer_name, pos_idx):
    """从稀疏或稠密矩阵中提取特定列的数据"""
    if layer_name not in adata.layers:
        # 如果碱基不是 ACGT，返回全 0 数组
        return np.zeros(adata.n_obs)
    
    data = adata.layers[layer_name][:, pos_idx]
    if issparse(data):
        return data.toarray().flatten()
    return data.flatten()

def main():
    parser = argparse.ArgumentParser(description="结合 SNV VarAllele 与 h5ad 全矩阵计算突变统计量")
    parser.add_argument("-s", "--snv", required=True, help="合并后的 SNV tsv 文件")
    parser.add_argument("-a", "--h5ad", required=True, help="全位点覆盖度的 h5ad 文件")
    parser.add_argument("-o", "--output", default="final_mutation_report.tsv", help="输出文件名")
    args = parser.parse_args()

    # 1. 加载数据
    print("正在加载数据...")
    # 自动识别分隔符（可能是空格或 Tab）
    snv_df = pd.read_csv(args.snv, sep='\t')
    adata = ad.read_h5ad(args.h5ad)
    
    var_to_idx = {name: i for i, name in enumerate(adata.var_names)}
    total_cells = adata.n_obs

    # 2. 筛选候选突变 (使用 VarAllele 作为真正的 Alt)
    print("正在从 VarAllele 筛选候选突变...")
    # 统计 (位置, 参考, 突变) 组合出现的细胞数
    mutation_groups = snv_df.groupby(['Position', 'Ref', 'VarAllele']).size().reset_index(name='cell_count')
    
    # 筛选在多于 1 个细胞中出现的突变
    candidate_muts = mutation_groups[mutation_groups['cell_count'] > 1].copy()
    print(f"筛选出 {len(candidate_muts)} 个候选突变位点 (Confident Cell > 1)。")

    results = []

    # 3. 遍历候选位点
    for _, row in candidate_muts.iterrows():
        pos = int(row['Position'])
        ref = str(row['Ref']).upper()
        alt = str(row['VarAllele']).upper() # 使用 VarAllele
        
        # 匹配 h5ad 坐标，兼容 "chrM:pos" 格式
        var_id = f"chrM:{pos}"
        if var_id not in var_to_idx:
            continue
        
        pos_idx = var_to_idx[var_id]

        # 提取全局数据 (2137个细胞)
        ref_fwd = get_layer_data(adata, f"{ref}_fwd", pos_idx)
        ref_rev = get_layer_data(adata, f"{ref}_rev", pos_idx)
        alt_fwd = get_layer_data(adata, f"{alt}_fwd", pos_idx)
        alt_rev = get_layer_data(adata, f"{alt}_rev", pos_idx)
        
        # 获取深度 X
        depth = adata.X[:, pos_idx].toarray().flatten() if issparse(adata.X) else adata.X[:, pos_idx]

        # --- 计算统计量 ---
        alt_total_per_cell = alt_fwd + alt_rev
        
        # VAF > 0 的细胞
        vaf_mask = (alt_total_per_cell > 0) & (depth > 0)
        vaf_pos_cells = np.sum(vaf_mask)
        
        if vaf_pos_cells > 0:
            # 只对有突变的细胞计算 VAF 统计量
            vaf_values = alt_total_per_cell[vaf_mask] / depth[vaf_mask]
            mean_vaf = np.mean(vaf_values)
            var_vaf = np.var(vaf_values)
            lis = mean_vaf / (1 + var_vaf)
        else:
            mean_vaf, var_vaf, lis = 0, 0, 0

        # Confident 定义: 突变碱基正负链各至少 1 条 read
        conf_mask = (alt_fwd >= 1) & (alt_rev >= 1)
        conf_cell_count = np.sum(conf_mask)

        alt_fw_total = alt_fwd.sum()
        alt_rev_total = alt_rev.sum()
        strand_score = alt_fw_total / (alt_fw_total + alt_rev_total) if (alt_fw_total + alt_rev_total) > 0 else 0

        results.append({
            'Position': pos,
            'Ref': ref,
            'Alt': alt,
            'Ref_fw_total': int(ref_fwd.sum()),
            'Ref_rev_total': int(ref_rev.sum()),
            'Alt_fw_total': int(alt_fw_total),
            'Alt_rev_total': int(alt_rev_total),
            'Strand_score': round(strand_score, 4),
            'Mean_vaf': round(mean_vaf, 6),
            'Var_vaf': round(var_vaf, 6),
            'Lis': round(lis, 6),
            'Pct_conf': round(conf_cell_count / vaf_pos_cells, 4),
            'Pct_vaf_pos': round(vaf_pos_cells / total_cells, 4),
            'N_cells_vaf_pos': int(vaf_pos_cells),
            'Conf_cells': int(conf_cell_count),
            'Total_cells': total_cells
        })

    # 4. 保存
    if not results:
        print("警告: 未生成任何统计结果，请检查 Chrom 命名是否匹配。")
        return

    res_df = pd.DataFrame(results)
    res_df = res_df.sort_values('Position')
    res_df.to_csv(args.output, sep='\t', index=False)
    print(f"完成！结果已保存至: {args.output}")

if __name__ == "__main__":
    main()