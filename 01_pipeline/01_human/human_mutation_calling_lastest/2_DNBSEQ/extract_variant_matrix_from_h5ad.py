import pandas as pd
import numpy as np
import anndata as ad
import gzip
import argparse
import os
from tqdm import tqdm
from scipy.sparse import issparse

def main():
    parser = argparse.ArgumentParser(description='从 h5ad 文件中提取指定 somatic 突变位点的计数矩阵')
    parser.add_argument('-s', '--somatic', required=True, help='筛选得到的 somatic 突变列表 (tsv)')
    parser.add_argument('-a', '--h5ad', required=True, help='All_cell_combined_counts_v2.h5ad 文件路径')
    parser.add_argument('-o', '--output', default='variant_sparse_matrix.tsv.gz', help='输出文件名')
    parser.add_argument('--chrom', default='chrM', help='染色体名称 (默认: chrM)')
    args = parser.parse_args()

    # 1. 加载突变列表
    print(f"正在读取突变列表: {args.somatic}")
    somatic_df = pd.read_csv(args.somatic, sep='\t')
    # 提取唯一的突变组合
    mut_list = somatic_df[['Position', 'Ref', 'Alt']].drop_duplicates().values.tolist()
    print(f"共需提取 {len(mut_list)} 个突变位点")

    # 2. 加载 AnnData
    print(f"正在加载 h5ad 文件: {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    
    # 构建 var 索引映射 (假设 var_names 格式为 "chrM:pos")
    var_names = adata.var_names.tolist()
    var_to_idx = {name: i for i, name in enumerate(var_names)}
    
    cells = adata.obs_names.tolist()
    n_cells = len(cells)

    # 3. 准备输出文件
    print(f"正在提取数据并写入 {args.output}...")
    with gzip.open(args.output, 'wt') as f:
        # 写入表头
        f.write("cell\tchrom\tpos\tref_base\talt_base\tref_count\talt_count\tvaf\n")
        
        # 遍历每一个突变
        for pos, ref, alt in tqdm(mut_list, desc="处理突变位点"):
            var_id = f"{args.chrom}:{int(pos)}"
            
            if var_id not in var_to_idx:
                # 如果位点不在 h5ad 中，跳过或记录
                continue
                
            idx = var_to_idx[var_id]
            
            # 获取 Ref 和 Alt 的计数 (从 layers 中提取)
            # 注意：这里合并了正负链
            ref_layer_f = f"{ref}_fwd"
            ref_layer_r = f"{ref}_rev"
            alt_layer_f = f"{alt}_fwd"
            alt_layer_r = f"{alt}_rev"
            
            # 提取所有细胞的该位点数据
            def get_data(layer_name):
                if layer_name not in adata.layers:
                    return np.zeros(n_cells)
                data = adata.layers[layer_name][:, idx]
                return data.toarray().flatten() if issparse(data) else data.flatten()

            ref_counts = get_data(ref_layer_f) + get_data(ref_layer_r)
            alt_counts = get_data(alt_layer_f) + get_data(alt_layer_r)
            
            # 获取深度 (X 矩阵通常存的是总 depth)
            depths = adata.X[:, idx]
            if issparse(depths):
                depths = depths.toarray().flatten()
            else:
                depths = depths.flatten()
                
            # 写入每一个细胞的记录
            for i in range(n_cells):
                r_c = int(ref_counts[i])
                a_c = int(alt_counts[i])
                total = depths[i]
                vaf = a_c / total if total > 0 else 0.0
                
                # 只有当有 reads 覆盖时才记录（为了真正的“稀疏”可以加这个条件，
                # 但根据你给的样例，全 0 似乎也要保留，这里按你样例全写出）
                line = f"{cells[i]}\t{args.chrom}\t{int(pos)}\t{ref}\t{alt}\t{r_c}\t{a_c}\t{vaf:.6f}\n"
                f.write(line)

    print(f"\n提取完成！生成文件大小: {os.path.getsize(args.output) / 1024 / 1024:.2f} MB")

if __name__ == "__main__":
    main()