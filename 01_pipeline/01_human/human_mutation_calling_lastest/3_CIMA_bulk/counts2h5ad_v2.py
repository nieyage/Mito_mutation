import pandas as pd
import numpy as np
import os
import glob
import argparse
from scipy.sparse import csr_matrix
import anndata as ad

def main():
    parser = argparse.ArgumentParser(description="将单细胞 .counts 文件合并为 H5AD 稀疏矩阵格式")
    parser.add_argument("-i", "--input", required=True, help="输入文件夹路径 (包含 .count 文件)")
    parser.add_argument("-o", "--output", required=True, help="输出 H5AD 文件路径 (例如 ./result.h5ad)")
    args = parser.parse_args()

    # 1. 获取所有文件
    file_list = sorted(glob.glob(os.path.join(args.input, "*.counts")))
    if not file_list:
        print(f"错误: 在 {args.input} 中未找到 .counts 文件")
        return

    print(f"找到 {len(file_list)} 个细胞文件，开始初始化...")

    # 2. 读取第一个文件以建立坐标参考 (Index)
    first_df = pd.read_csv(file_list[0], sep='\t')
    # 使用 chrom_pos 作为唯一位点标识
    var_names = (first_df['chrom'] + ":" + first_df['pos'].astype(str)).values
    var_df = first_df[['chrom', 'pos', 'ref']].copy()
    var_df.index = var_names
    
    n_vars = len(var_names)
    n_cells = len(file_list)
    cell_names = []

    # 3. 准备存储稀疏矩阵数据的字典
    # 我们将 depth 存为 X，其他碱基计数存为 layers
    columns_to_extract = ['depth', 'A_fwd', 'A_rev', 'C_fwd', 'C_rev', 'G_fwd', 'G_rev', 'T_fwd', 'T_rev']
    data_dict = {col: [] for col in columns_to_extract}
    row_indices = {col: [] for col in columns_to_extract}
    col_indices = {col: [] for col in columns_to_extract}

    # 4. 遍历文件读取数据
    for cell_idx, file_path in enumerate(file_list):
        cell_name = os.path.basename(file_path).replace(".counts", "")
        cell_names.append(cell_name)
        
        # 读取当前细胞数据
        df = pd.read_csv(file_path, sep='\t')
        
        for col in columns_to_extract:
            # 只提取非零值以构建稀疏矩阵
            mask = df[col] > 0
            if mask.any():
                vals = df.loc[mask, col].values
                idxs = np.where(mask)[0]
                
                data_dict[col].extend(vals)
                row_indices[col].extend([cell_idx] * len(vals)) # AnnData 默认 (cell, var)
                col_indices[col].extend(idxs)

        if (cell_idx + 1) % 100 == 0:
            print(f"进度: {cell_idx + 1}/{n_cells} 细胞已加载...")

    print("正在构建 AnnData 对象...")

    # 5. 构建稀疏矩阵并组装 AnnData
    # 注意：AnnData 的标准形状是 (n_obs, n_vars)，即 (细胞, 位点)
    obs_df = pd.DataFrame(index=cell_names)
    
    # 将 depth 作为主矩阵 X
    X_matrix = csr_matrix(
        (data_dict['depth'], (row_indices['depth'], col_indices['depth'])),
        shape=(n_cells, n_vars)
    )

    adata = ad.AnnData(X=X_matrix, obs=obs_df, var=var_df)

    # 将其他碱基计数作为 layers 存储
    for col in columns_to_extract:
        if col == 'depth': continue
        adata.layers[col] = csr_matrix(
            (data_dict[col], (row_indices[col], col_indices[col])),
            shape=(n_cells, n_vars)
        )

    # 6. 保存文件
    print(f"正在写入文件: {args.output} ...")
    adata.write_h5ad(args.output)
    print("转换成功！")

if __name__ == "__main__":
    main()