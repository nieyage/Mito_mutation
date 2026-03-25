import pandas as pd
import numpy as np
import anndata as ad
import argparse
from scipy.sparse import issparse
import sys

def main():
    parser = argparse.ArgumentParser(description="高性能线粒体突变统计量计算 (支持大文件 & 自动过滤 INDEL)")
    parser.add_argument("-s", "--snv", required=True, help="合并后的 SNV tsv 文件")
    parser.add_argument("-a", "--h5ad", required=True, help="全位点覆盖度的 h5ad 文件")
    parser.add_argument("-o", "--output", default="final_mutation_report.tsv", help="输出文件名")
    args = parser.parse_args()

    # 1. 加载并清洗 SNV 数据
    print("正在加载并清洗 SNV 数据...")
    snv_df = pd.read_csv(args.snv, sep='\t')
    initial_count = len(snv_df)
    
    # 关键步骤：只保留 A, C, G, T 的标准替换 (SNV)
    valid_bases = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}
    snv_df = snv_df[
        snv_df['Ref'].isin(valid_bases) & 
        snv_df['VarAllele'].isin(valid_bases)
    ].copy()
    snv_df['Ref'] = snv_df['Ref'].str.upper()
    snv_df['VarAllele'] = snv_df['VarAllele'].str.upper()
    
    print(f"原始记录: {initial_count}, 过滤 INDEL 后剩余 SNV: {len(snv_df)}")

    # 2. 初始化 h5ad (使用 backed 模式节省内存)
    print("正在初始化 h5ad 文件索引...")
    adata = ad.read_h5ad(args.h5ad, backed='r')
    total_cells = adata.n_obs
    var_to_idx = {name: i for i, name in enumerate(adata.var_names)}

    # 筛选在多于 1 个细胞中出现的突变
    mutation_groups = snv_df.groupby(['Position', 'Ref', 'VarAllele']).size().reset_index(name='cell_count')
    candidate_muts = mutation_groups[mutation_groups['cell_count'] > 1].copy()
    
    # 映射位点索引
    candidate_muts['pos_idx'] = candidate_muts['Position'].apply(lambda x: var_to_idx.get(f"chrM:{x}", -1))
    candidate_muts = candidate_muts[candidate_muts['pos_idx'] != -1]
    
    unique_pos_indices = sorted(candidate_muts['pos_idx'].unique())
    print(f"待处理有效位点数: {len(unique_pos_indices)}")

    if not unique_pos_indices:
        print("错误: 未找到匹配的位点，请检查 Position 格式。")
        return

    # 3. 批量提取数据 (最耗时但最高效)
    print("正在批量提取相关位点数据至内存 (这可能需要几分钟)...")
    
    def fetch_matrix_cols(layer_key, indices):
        # 内部函数：高效提取特定列
        try:
            if layer_key is None:
                mat = adata.X[:, indices]
            elif layer_key in adata.layers:
                mat = adata.layers[layer_key][:, indices]
            else:
                return np.zeros((total_cells, len(indices)))
            
            # 如果是存储在磁盘的稀疏矩阵，toarray() 会将其载入内存
            return mat.toarray() if issparse(mat) else np.array(mat)
        except Exception as e:
            print(f"提取层 {layer_key} 时出错: {e}")
            return np.zeros((total_cells, len(indices)))

    # 加载 A, C, G, T 的 fwd/rev 和 Depth
    loaded_data = {}
    idx_map = {orig_idx: i for i, orig_idx in enumerate(unique_pos_indices)}
    
    bases = ['A', 'C', 'G', 'T']
    for b in bases:
        for suffix in ['_fwd', '_rev']:
            layer_name = f"{b}{suffix}"
            loaded_data[layer_name] = fetch_matrix_cols(layer_name, unique_pos_indices)
    
    depth_matrix = fetch_matrix_cols(None, unique_pos_indices)

    # 4. 向量化统计计算
    print("正在进行向量化统计计算...")
    results = []
    
    for _, row in candidate_muts.iterrows():
        pos = int(row['Position'])
        ref = row['Ref']
        alt = row['VarAllele']
        col_i = idx_map[row['pos_idx']]

        # 从预加载的内存数据中切片
        rf = loaded_data[f"{ref}_fwd"][:, col_i]
        rr = loaded_data[f"{ref}_rev"][:, col_i]
        af = loaded_data[f"{alt}_fwd"][:, col_i]
        ar = loaded_data[f"{alt}_rev"][:, col_i]
        d  = depth_matrix[:, col_i]

        alt_total_per_cell = af + ar
        
        # VAF 统计
        vaf_mask = (alt_total_per_cell > 0) & (d > 0)
        vaf_pos_cells = np.sum(vaf_mask)
        
        if vaf_pos_cells > 0:
            vaf_values = alt_total_per_cell[vaf_mask] / d[vaf_mask]
            mean_vaf = np.mean(vaf_values)
            var_vaf = np.var(vaf_values)
            lis = mean_vaf / (1 + var_vaf)
        else:
            mean_vaf = var_vaf = lis = 0

        # Confident 定义 (正负链各 >= 1)
        conf_mask = (af >= 1) & (ar >= 1)
        conf_cell_count = np.sum(conf_mask)

        af_sum, ar_sum = af.sum(), ar.sum()
        strand_score = af_sum / (af_sum + ar_sum) if (af_sum + ar_sum) > 0 else 0

        results.append({
            'Position': pos,
            'Ref': ref,
            'Alt': alt,
            'Ref_fw_total': int(rf.sum()),
            'Ref_rev_total': int(rr.sum()),
            'Alt_fw_total': int(af_sum),
            'Alt_rev_total': int(ar_sum),
            'Strand_score': round(strand_score, 4),
            'Mean_vaf': round(mean_vaf, 6),
            'Var_vaf': round(var_vaf, 6),
            'Lis': round(lis, 6),
            'Pct_conf': round(conf_cell_count / vaf_pos_cells, 4) if vaf_pos_cells > 0 else 0,
            'Pct_vaf_pos': round(vaf_pos_cells / total_cells, 4),
            'N_cells_vaf_pos': int(vaf_pos_cells),
            'Conf_cells': int(conf_cell_count),
            'Total_cells': total_cells
        })

    # 5. 保存结果
    if not results:
        print("警告: 无结果生成。")
        return

    res_df = pd.DataFrame(results).sort_values('Position')
    res_df.to_csv(args.output, sep='\t', index=False)
    print(f"完成！结果已保存至: {args.output}")

if __name__ == "__main__":
    main()