import pandas as pd
import argparse
import sys

def run_mitochondrial_filter(input_file, output_file, min_vaf, min_reads):
    try:
        # 1. 读取数据
        # 自动识别分隔符（针对空格或 Tab 对齐的情况）
        df = pd.read_csv(input_file, sep='\s+')
    except Exception as e:
        print(f"读取文件失败: {e}")
        sys.exit(1)

    print(f"开始处理: {input_file}")
    original_count = len(df)
    
    # --- 格式处理 ---
    # 处理 VarFreq：去掉 '%' 并转为小数 (例如 "0.28%" -> 0.0028)
    if df['VarFreq'].dtype == 'O':
        df['VarFreq'] = df['VarFreq'].str.rstrip('%').astype('float') / 100.0

    # --- 过滤条件 1: 基础频率 ---
    df = df[df['VarFreq'] >= min_vaf]

    # --- 过滤条件 2: 排除特定错配区域 (Misalignment) ---
    # 定义排除列表 (Position, VarAllele)
    exclude_list = [
        (302, 'C'), (309, 'T'), (311, 'T'), (312, 'T'), (313, 'T'), (316, 'C'),
        (514, 'A'), (515, 'G'), (523, 'C'), (524, 'G'), (3107, 'ANY')
    ]
    for pos, allele in exclude_list:
        if allele == 'ANY':
            df = df[df['Position'] != pos]
        else:
            df = df[~((df['Position'] == pos) & (df['VarAllele'] == allele))]

    # --- 过滤条件 3: 支撑 Reads 数与链对称性 ---
    # 要求：正向和反向 Reads 均需 > 指定值 (默认 2)
    df = df[(df['Reads2Plus'] > min_reads) & (df['Reads2Minus'] > min_reads)]
    
    # 比例：30% < forward / (forward + reverse) < 70%
    df['Strand_Ratio'] = df['Reads2Plus'] / (df['Reads2Plus'] + df['Reads2Minus'])
    df = df[(df['Strand_Ratio'] > 0.3) & (df['Strand_Ratio'] < 0.7)]
    df = df.rename(columns={
        'cell_id': 'celltype',
        'VarAllele': 'Alt'
    })

    # --- 过滤条件 4: 群体高频突变 (0.9 < VAF < 1.0) ---
    df = df[df['VarFreq'] < 0.9]

    # 保存结果
    df.to_csv(output_file, sep='\t', index=False)
    print(f"过滤完成！")
    print(f"原始条数: {original_count}")
    print(f"保留条数: {len(df)}")
    print(f"结果已保存至: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="线粒体变异筛选脚本 (Cell-type Level)")
    
    # 添加参数
    parser.add_argument("-i", "--input", required=True, help="输入的 SNV 表格文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出的过滤后 CSV 文件路径")
    parser.add_argument("-v", "--vaf", type=float, default=0.001, help="最小变异频率阈值 (默认: 0.001)")
    parser.add_argument("-r", "--reads", type=int, default=2, help="最小单向支撑 reads 数 (默认: >2)")

    args = parser.parse_args()

    # 运行过滤函数
    run_mitochondrial_filter(args.input, args.output, args.vaf, args.reads)