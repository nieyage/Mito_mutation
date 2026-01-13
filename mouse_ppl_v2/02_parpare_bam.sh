# 批量处理Step1-5：
cut -d, -f1 /md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv | tr -d '"' > barcodes_BMMC27m.txt
cut -d" " -f1 /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/somatic_mutation/AR3_C4_celltype_info.txt | tr -d '"' > barcodes_AR3_C4.txt
cut -d" " -f1 /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/somatic_mutation/AR3_C5_celltype_info.txt | tr -d '"' > barcodes_AR3_C5.txt


# 1. sample_list.txt
# 格式：样本名 BAM文件路径 Barcode文件路径

BMMC_27m /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/BMMC_27m_ATAC/outs/atac_possorted_bam.bam barcodes_BMMC27m.txt
AR3_C4 /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/outs/possorted_bam.bam barcodes_AR3_C4.txt
AR3_C5 /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/outs/possorted_bam.bam barcodes_AR3_C5.txt
# 2. batch_process_mito.sh

#!/bin/bash
# 批量处理脚本：batch_process_mito.sh
# 用法：./batch_process_mito.sh <样本列表文件>

SAMPLE_LIST="$1"

if [ -z "$SAMPLE_LIST" ]; then
    echo "用法: $0 <样本列表文件>"
    echo "样本列表文件格式："
    echo "样本名1 /path/to/bam1.bam /path/to/barcodes1.txt"
    echo "样本名2 /path/to/bam2.bam /path/to/barcodes2.txt"
    exit 1
fi

while IFS= read -r line || [[ -n "$line" ]]; do
    # 跳过空行和注释
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    
    # 解析行
    read -r sample_name bam_file barcode_file <<< "$line"
    
    echo "处理样本: $sample_name"
    echo "BAM文件: $bam_file"
    echo "Barcode文件: $barcode_file"
    
    # 运行处理脚本
    ./process_mito_bam.sh "$bam_file" "$sample_name" "$barcode_file" 36
    
    echo "样本 $sample_name 处理完成"
    echo "---------------------------"
done < "$SAMPLE_LIST"

echo "所有样本处理完成！"


# 3. process_mito_bam.sh

#!/bin/bash
# 脚本名称：process_mito_bam.sh
# 用途：处理cellranger输出的BAM文件，提取线粒体相关reads并重新比对到偏移的chrM基因组
# 用法：./process_mito_bam.sh <input_bam> <output_prefix> <barcode_file> [threads]

set -e  # 遇到错误时退出

# 参数检查
if [ $# -lt 3 ]; then
    echo "用法: $0 <input_bam> <output_prefix> <barcode_file> [threads]"
    echo "示例: $0 /path/to/atac_possorted_bam.bam PBMC /path/to/barcodes.txt 36"
    exit 1
fi

# 参数设置
INPUT_BAM="$1"
OUTPUT_PREFIX="$2"
BARCODE_FILE="$3"
THREADS="${4:-36}"  # 默认36线程

# 路径设置（根据实际情况修改）
SHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta"
UNSHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"

echo "=========================================="
echo "开始处理: $INPUT_BAM"
echo "输出前缀: $OUTPUT_PREFIX"
echo "线程数: $THREADS"
echo "=========================================="

# 创建输出目录
mkdir -p ${OUTPUT_PREFIX}_output/unshifted_bam
mkdir -p ${OUTPUT_PREFIX}_output/Dloop_bam
mkdir -p ${OUTPUT_PREFIX}_output/tmp

# Step 1: 提取chrM reads和未比对reads
echo "Step 1: 提取chrM reads和未比对reads..."
samtools view -@ $THREADS -b $INPUT_BAM chrM > ${OUTPUT_PREFIX}_output/chrM_mapped.bam
samtools view -@ $THREADS -b -f 4 $INPUT_BAM > ${OUTPUT_PREFIX}_output/unmapped.bam
samtools merge -@ $THREADS ${OUTPUT_PREFIX}_output/possorted_chrM_and_unmapped.bam \
    ${OUTPUT_PREFIX}_output/unmapped.bam ${OUTPUT_PREFIX}_output/chrM_mapped.bam
samtools index -@ $THREADS ${OUTPUT_PREFIX}_output/possorted_chrM_and_unmapped.bam
rm ${OUTPUT_PREFIX}_output/unmapped.bam

# Step 2: 使用偏移的chrM基因组重新比对
echo "Step 2: 使用偏移的chrM基因组重新比对..."
PBMC_CHRM_UNMAPPED_BAM="${OUTPUT_PREFIX}_output/possorted_chrM_and_unmapped.bam"

samtools collate -@ $THREADS -Oun128 $PBMC_CHRM_UNMAPPED_BAM \
    | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
    | bwa mem -pt$THREADS -CH <(samtools view -H $PBMC_CHRM_UNMAPPED_BAM | grep ^@RG) $SHIFTED_CHRM_REF - \
    | samtools sort -@$THREADS -m4g -o ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted.bam -

# 过滤未比对reads
echo "Step 2.1: 过滤和索引..."
samtools index -@ $THREADS ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted.bam
samtools view -@ $THREADS -b ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted.bam chrM > ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted2.bam

# 只保留D-loop区域 (chrM:7000-10000)
samtools index -@ $THREADS ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted2.bam
samtools view -@ $THREADS -h -b ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted2.bam chrM:7000-10000 > ${OUTPUT_PREFIX}_output/chrM_Dloop.bam

# Step 3: 按CB标签排序
echo "Step 3: 按CB标签排序..."
ALL_BAM="${OUTPUT_PREFIX}_output/chrM_mapped.bam"
DLOOP_BAM="${OUTPUT_PREFIX}_output/chrM_Dloop.bam"

samtools sort -@ $THREADS -t CB $ALL_BAM -o ${OUTPUT_PREFIX}_output/chrM_sorted_CB.bam
samtools sort -@ $THREADS -t CB $DLOOP_BAM -o ${OUTPUT_PREFIX}_output/chrM_Dloop_sorted_CB.bam

# Step 4: 根据barcode拆分BAM文件
echo "Step 4: 根据barcode拆分BAM文件..."

# 运行拆分脚本
echo "拆分未偏移的BAM文件..."
python split_bam.py -i ${OUTPUT_PREFIX}_output/chrM_sorted_CB.bam \
    -b $BARCODE_FILE \
    -o ${OUTPUT_PREFIX}_output/unshifted_bam/

echo "拆分D-loop BAM文件..."
python split_bam.py -i ${OUTPUT_PREFIX}_output/chrM_Dloop_sorted_CB.bam \
    -b $BARCODE_FILE \
    -o ${OUTPUT_PREFIX}_output/Dloop_bam/

# 索引拆分后的BAM文件（可选）
echo "索引拆分后的BAM文件..."
for bam in ${OUTPUT_PREFIX}_output/unshifted_bam/*.bam ${OUTPUT_PREFIX}_output/Dloop_bam/*.bam; do
    samtools index -@ $THREADS $bam &
done
wait

# 清理中间文件（可选）
echo "清理中间文件..."
rm -f ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted.bam.bai
rm -f ${OUTPUT_PREFIX}_output/chrM_unmapped_bwa_shifted2.bam.bai

echo "=========================================="
echo "处理完成！"
echo "输出目录: ${OUTPUT_PREFIX}_output/"
echo "未偏移BAM: ${OUTPUT_PREFIX}_output/unshifted_bam/"
echo "D-loop BAM: ${OUTPUT_PREFIX}_output/Dloop_bam/"
echo "=========================================="

chmod +x batch_process_mito.sh process_mito_bam.sh

nohup ./batch_process_mito.sh sample_list.txt > bam_processing.log 2>&1 &

# Step6: SNV calling (注意！！这一步不要多个样本一起跑！！！)

./process_mito_variants.sh <csv_file> <output_base> <unshifted_bam_base> <shifted_bam_base>

./process_mito_variants.sh \
    /path/to/barcodes.csv \
    /path/to/output \
    /path/to/unshifted_bams \
    /path/to/shifted_bams


./process_mito_variants.sh \
    barcodes_BMMC27m.txt \
    ./BMMC_27m_output/SNVcalling_allcell \
    ./BMMC_27m_output/unshifted_bam/ \
    ./BMMC_27m_output/Dloop_bam/

./process_mito_variants.sh \
    barcodes_AR3_C4.txt \
    ./AR3_C4_output/SNVcalling_allcell \
    ./AR3_C4_output/unshifted_bam/ \
    ./AR3_C4_output/Dloop_bam/

./process_mito_variants.sh \
    barcodes_AR3_C5.txt \
    ./AR3_C5_output/SNVcalling_allcell \
    ./AR3_C5_output/unshifted_bam/ \
    ./AR3_C5_output/Dloop_bam/

nohup sh ./snvcalling_run.sh > SNV_calling_running.log 2>&1 &
