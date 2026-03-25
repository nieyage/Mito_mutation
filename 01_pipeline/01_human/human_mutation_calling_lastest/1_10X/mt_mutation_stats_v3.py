import pandas as pd
import numpy as np
import anndata as ad
import argparse
from scipy.sparse import issparse
import sys

def main():
    parser = argparse.ArgumentParser(description="线粒体突变统计 (Confident 判断基于 TSV 中的 Reads2Plus/Minus)")
    parser.add_argument("-s", "--snv", required=True, help="合并后的 SNV tsv 文件")
    parser.add_argument("-a", "--h5ad", required=True, help="全位点覆盖度的 h5ad 文件")
    parser.add_argument("-o", "--output", default="final_mutation_report.tsv", help="输出文件名")
    args = parser.parse_args()

    # 1. 加载并清洗 SNV 数据
    print("正在加载并清洗 SNV 数据...")
    snv_df = pd.read_csv(args.snv, sep='\t')
    
    # 检查必要列
    required_cols = ['VarFreq', 'Reads2Plus', 'Reads2Minus', 'Ref', 'VarAllele', 'Position']
    for col in required_cols:
        if col not in snv_df.columns:
            print(f"错误：在 SNV 文件中未找到 '{col}' 列！")
            sys.exit(1)

    # A. 处理 VarFreq (3.03% -> 0.0303)
    snv_df['VarFreq'] = snv_df['VarFreq'].astype(str).str.replace('%', '', regex=False).astype(float) / 100.0
    
    # B. 关键：在行级别判断该细胞是否是 Confident
    # 标准：突变碱基的正链和负链 read 数均 >= 1
    snv_df['is_conf'] = (snv_df['Reads2Plus'] >= 1) & (snv_df['Reads2Minus'] >= 1)
    
    # C. 过滤 INDEL，保留标准 SNV
    valid_bases = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}
    snv_df = snv_df[
        snv_df['Ref'].isin(valid_bases) & 
        snv_df['VarAllele'].isin(valid_bases)
    ].copy()
    snv_df['Ref'] = snv_df['Ref'].str.upper()
    snv_df['VarAllele'] = snv_df['VarAllele'].str.upper()

    # 2. 初始化 h5ad
    print("正在初始化 h5ad 文件索引...")
    adata = ad.read_h5ad(args.h5ad, backed='r')
    total_cells = adata.n_obs
    var_to_idx = {name: i for i, name in enumerate(adata.var_names)}

    # --- 核心修改：聚合统计 ---
    # 同时计算：总细胞数、Confident细胞数、VAF均值、VAF方差
    mutation_groups = snv_df.groupby(['Position', 'Ref', 'VarAllele']).agg(
        total_recorded_cells=('Position', 'size'),           # TSV 中记录的所有细胞
        conf_cell_count=('is_conf', 'sum'),                  # 满足正负链支持的细胞
        mean_vaf_pre=('VarFreq', 'mean'),
        var_vaf_pre=('VarFreq', lambda x: np.var(x, ddof=0))
    ).reset_index()
    
    # 筛选在多于 1 个细胞中出现的突变
    candidate_muts = mutation_groups[mutation_groups['total_recorded_cells'] > 1].copy()
    
    # 映射位点索引
    candidate_muts['pos_idx'] = candidate_muts['Position'].apply(lambda x: var_to_idx.get(f"chrM:{x}", -1))
    candidate_muts = candidate_muts[candidate_muts['pos_idx'] != -1]
    
    unique_pos_indices = sorted(candidate_muts['pos_idx'].unique())
    print(f"待处理有效位点数: {len(unique_pos_indices)}")

    # 3. 批量提取数据 (仅用于计算 Ref_total 和 Strand_score 等背景指标)
    print("正在批量提取 h5ad 全局数据...")
    idx_map = {orig_idx: i for i, orig_idx in enumerate(unique_pos_indices)}
    
    def fetch_matrix_cols(layer_key, indices):
        if layer_key is None: mat = adata.X[:, indices]
        elif layer_key in adata.layers: mat = adata.layers[layer_key][:, indices]
        else: return np.zeros((total_cells, len(indices)))
        return mat.toarray() if issparse(mat) else np.array(mat)

    loaded_data = {}
    for b in ['A', 'C', 'G', 'T']:
        for suffix in ['_fwd', '_rev']:
            l = f"{b}{suffix}"
            loaded_data[l] = fetch_matrix_cols(l, unique_pos_indices)
    depth_matrix = fetch_matrix_cols(None, unique_pos_indices)

    # 4. 统计计算
    print("正在计算最终统计量...")
    results = []
    for _, row in candidate_muts.iterrows():
        pos, ref, alt = int(row['Position']), row['Ref'], row['VarAllele']
        col_i = idx_map[row['pos_idx']]

        # 获取全局 read 计数（用于 Strand_score）
        af = loaded_data[f"{alt}_fwd"][:, col_i]
        ar = loaded_data[f"{alt}_rev"][:, col_i]
        rf = loaded_data[f"{ref}_fwd"][:, col_i]
        rr = loaded_data[f"{ref}_rev"][:, col_i]
        d  = depth_matrix[:, col_i]

        # 计算 VAF 检测背景
        vaf_mask = ((af + ar) > 0) & (d > 0)
        vaf_pos_cells = np.sum(vaf_mask)
        
        # 提取聚合后的统计值
        conf_cells = int(row['conf_cell_count'])
        mean_vaf = row['mean_vaf_pre'] if vaf_pos_cells > 0 else 0
        var_vaf = row['var_vaf_pre'] if vaf_pos_cells > 0 else 0
        lis = mean_vaf / (1 + var_vaf) if vaf_pos_cells > 0 else 0

        af_sum, ar_sum = af.sum(), ar.sum()

        results.append({
            'Position': pos,
            'Ref': ref,
            'Alt': alt,
            'Ref_fw_total': int(rf.sum()),
            'Ref_rev_total': int(rr.sum()),
            'Alt_fw_total': int(af_sum),
            'Alt_rev_total': int(ar_sum),
            'Strand_score': round(af_sum / (af_sum + ar_sum), 4) if (af_sum + ar_sum) > 0 else 0,
            'Mean_vaf': round(mean_vaf, 6),
            'Var_vaf': round(var_vaf, 6),
            'Lis': round(lis, 6),
            'Pct_conf': round(conf_cells / vaf_pos_cells, 4) if vaf_pos_cells > 0 else 0,
            'Pct_vaf_pos': round(vaf_pos_cells / total_cells, 4),
            'N_cells_vaf_pos': int(vaf_pos_cells),
            'Conf_cells': conf_cells,
            'Total_cells': total_cells
        })

    # 5. 保存
    pd.DataFrame(results).sort_values('Position').to_csv(args.output, sep='\t', index=False)
    print(f"完成！结果保存至: {args.output}")

if __name__ == "__main__":
    main()