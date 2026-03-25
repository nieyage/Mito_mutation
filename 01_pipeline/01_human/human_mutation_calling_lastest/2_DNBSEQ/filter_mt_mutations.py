import pandas as pd
import numpy as np
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='根据统计指标筛选线粒体 Germline 和 Somatic 突变')
    
    # 输入输出
    parser.add_argument('-i', '--input', required=True, help='输入的汇总统计文件 (tsv)')
    parser.add_argument('-o', '--output-prefix', default='filtered', help='输出文件前缀')
    
    # 筛选参数
    parser.add_argument('--strand-min', type=float, default=0.2, help='Strand_score 最小值 (默认: 0.2)')
    parser.add_argument('--strand-max', type=float, default=0.8, help='Strand_score 最大值 (默认: 0.8)')
    parser.add_argument('--min-conf-cells', type=int, default=2, help='最小 Conf_cells 数量 (默认: 2)')
    parser.add_argument('--germline-vaf', type=float, default=0.9, help='判定为 Germline 的 Mean_vaf 阈值 (默认: 0.9)')
    parser.add_argument('--somatic-vaf-max', type=float, default=0.5, help='判定为 Somatic 的 Mean_vaf 上限 (默认: 0.5)')

    args = parser.parse_args()

    # 1. 读取数据
    print(f"读取文件: {args.input}")
    df = pd.read_csv(args.input, sep='\t')
    initial_count = len(df)

    # 2. 基础质量过滤 (所有突变都必须满足)
    # - 链偏好性过滤：双链都要有支持，防止测序伪影
    # - 细胞数量过滤：至少在指定数量的细胞中是 "Confident" 的
    quality_mask = (
        (df['Strand_score'] >= args.strand_min) & 
        (df['Strand_score'] <= args.strand_max) &
        (df['Conf_cells'] >= args.min_conf_cells)
    )
    
    df_q = df[quality_mask].copy()
    print(f"经过质量过滤 (Strand & Conf_cells >= {args.min_conf_cells})，剩余: {len(df_q)}/{initial_count}")

    # 3. 筛选 Germline 突变
    # 特征：均值 VAF 极高，几乎所有细胞都携带
    germline_df = df_q[df_q['Mean_vaf'] >= args.germline_vaf].copy()

    # 4. 筛选 Somatic 突变
    # 特征：均值 VAF 较低（异质性），且 Lis 得分通常较高（表示在阳性细胞中稳定）
    somatic_df = df_q[df_q['Mean_vaf'] <= args.somatic_vaf_max].copy()

    # 5. 保存结果
    germline_file = f"{args.output_prefix}_germline.tsv"
    somatic_file = f"{args.output_prefix}_somatic.tsv"
    
    germline_df.sort_values('Position').to_csv(germline_file, sep='\t', index=False)
    somatic_df.sort_values('Position').to_csv(somatic_file, sep='\t', index=False)

    # 6. 打印总结
    print("\n" + "="*30)
    print(f"筛选总结报告:")
    print(f"- 总输入突变位点: {initial_count}")
    print(f"- Germline 突变 (Mean_VAF >= {args.germline_vaf}): {len(germline_df)}")
    print(f"- Somatic 突变 (Mean_VAF <= {args.somatic_vaf_max}): {len(somatic_df)}")
    print(f"- 结果已保存至: {germline_file} 和 {somatic_file}")
    print("="*30)

    if not somatic_df.empty:
        print("\n高可信 Somatic 突变预览 (按 Lis 得分排序):")
        print(somatic_df.sort_values('Lis', ascending=False).head(10)[['Position', 'Ref', 'Alt', 'Mean_vaf', 'Lis', 'Conf_cells']])

if __name__ == "__main__":
    main()