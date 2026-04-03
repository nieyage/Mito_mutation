# for batch1 

# generate_list.py

import pandas as pd
import os

# 配置信息
SAMPLE_LIST_FILE = "/data/R01/renwx5/MT_mution/CIMA/outs/mt_ratio/MT_ratio_gt_0.1_1.csv"
METADATA_FILE = "/md01/renwx5/MT_mution/CIMA/metadata_CNP0006227_merged.csv"
BASE_DATA_DIR = "/data/backup/renwx5/CIMA_data/CNP0006227_1"
OUTPUT_FILE = "batch1_task_list.txt"

def main():
    # 读取目标样本名
    samples_df = pd.read_csv(SAMPLE_LIST_FILE)
    target_samples = samples_df['sample'].tolist()

    # 读取元数据
    meta_df = pd.read_csv(METADATA_FILE)
    
    tasks = []
    print(f"正在匹配路径...")

    for sample in target_samples:
        match = meta_df[meta_df['run_accession'] == sample]
        if not match.empty:
            s_acc = match.iloc[0]['sample_accession'] 
            e_acc = match.iloc[0]['experiment_accession']
            fastq_dir = os.path.join(BASE_DATA_DIR, s_acc, e_acc, sample)
            
            if os.path.exists(fastq_dir):
                tasks.append(f"{sample} {fastq_dir}")
            else:
                print(f"路径不存在: {fastq_dir}")
        else:
            print(f"元数据未找到样本: {sample}")

    # 写入列表文件
    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(tasks))
    
    print(f"成功生成列表: {OUTPUT_FILE}, 共 {len(tasks)} 个样本")

if __name__ == "__main__":
    main()


# run_dnbc4_batch1.pbs

#!/bin/bash
#PBS -N ATAC_batch1
#PBS -q batch
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -t 1-26%4

# 1. 进入工作目录
cd /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/batch_1

# 2. 定义变量
DNBC4TOOLS="/md01/nieyg/software/dnbc4tools3.0/dnbc4tools"
GENOME_DIR="/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_masked_ATAC_addflank"
LIST="batch1_task_list.txt"

# 3. 根据当前的 Array ID 获取对应的行内容
# PBS_ARRAYID 会从 1 自动变到 26
LINE=$(sed -n "${PBS_ARRAYID}p" $LIST)
SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FASTQ_PATH=$(echo $LINE | awk '{print $2}')

# 4. 打印开始信息到主日志
echo "Task $PBS_ARRAYID started for $SAMPLE_ID at $(date)"

# 5. 运行命令
$DNBC4TOOLS atac run \
  --name $SAMPLE_ID \
  --fastqs $FASTQ_PATH \
  --genomeDir $GENOME_DIR \
  --threads 16 \
  --need_bam

echo "Task $PBS_ARRAYID for $SAMPLE_ID finished at $(date)"







# generate_list_batch2.py

import pandas as pd
import os

# 配置信息
SAMPLE_LIST_FILE = "/data/R01/renwx5/MT_mution/CIMA/outs/mt_ratio/MT_ratio_gt_0.1_2.csv"
METADATA_FILE = "/md01/renwx5/MT_mution/CIMA/metadata_CNP0006227_merged.csv"
BASE_DATA_DIR = "/data/backup/renwx5/CIMA_data/CNP0006227_2"
OUTPUT_FILE = "batch2_task_list.txt"

def main():
    # 读取目标样本名
    samples_df = pd.read_csv(SAMPLE_LIST_FILE)
    target_samples = samples_df['sample'].tolist()

    # 读取元数据
    meta_df = pd.read_csv(METADATA_FILE)
    
    tasks = []
    print(f"正在匹配路径...")

    for sample in target_samples:
        match = meta_df[meta_df['run_accession'] == sample]
        if not match.empty:
            s_acc = match.iloc[0]['sample_accession'] 
            e_acc = match.iloc[0]['experiment_accession']
            fastq_dir = os.path.join(BASE_DATA_DIR, s_acc, e_acc, sample)
            
            if os.path.exists(fastq_dir):
                tasks.append(f"{sample} {fastq_dir}")
            else:
                print(f"路径不存在: {fastq_dir}")
        else:
            print(f"元数据未找到样本: {sample}")

    # 写入列表文件
    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(tasks))
    
    print(f"成功生成列表: {OUTPUT_FILE}, 共 {len(tasks)} 个样本")

if __name__ == "__main__":
    main()


# run_dnbc4_batch2.pbs

#!/bin/bash
#PBS -N ATAC_batch2
#PBS -q batch
#PBS -l nodes=1:ppn=64
#PBS -j oe
#PBS -t 1-19%4

# 1. 进入工作目录
cd /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/batch_2

# 2. 定义变量
DNBC4TOOLS="/md01/nieyg/software/dnbc4tools3.0/dnbc4tools"
GENOME_DIR="/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_masked_ATAC_addflank"
LIST="batch2_task_list.txt"

# 3. 根据当前的 Array ID 获取对应的行内容
# PBS_ARRAYID 会从 1 自动变到 26
LINE=$(sed -n "${PBS_ARRAYID}p" $LIST)
SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FASTQ_PATH=$(echo $LINE | awk '{print $2}')

# 4. 打印开始信息到主日志
echo "Task $PBS_ARRAYID started for $SAMPLE_ID at $(date)"

# 5. 运行命令
$DNBC4TOOLS atac run \
  --name $SAMPLE_ID \
  --fastqs $FASTQ_PATH \
  --genomeDir $GENOME_DIR \
  --threads 16 \
  --need_bam

echo "Task $PBS_ARRAYID for $SAMPLE_ID finished at $(date)"


