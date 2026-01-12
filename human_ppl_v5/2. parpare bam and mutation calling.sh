# ref:https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
# https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references
bwa index Homo_sapiens_assembly38.chrM.fasta
bwa index Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
samtools faidx Homo_sapiens_assembly38.chrM.fasta
samtools faidx Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta


# Step1: Get unmapped and chrM reads
# input_bam= /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
input_bam=/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam
samtools view -b $input_bam chrM > PBMC_chrM_mapped.bam
samtools view -b -f 4 $input_bam > PBMC_unmapped.bam
samtools merge PBMC_possorted_chrM_and_unmapped.bam PBMC_unmapped.bam PBMC_chrM_mapped.bam
samtools index PBMC_possorted_chrM_and_unmapped.bam
rm PBMC_unmapped.bam 

# Step2: Remapped with shifted chrM genome by bwa

#  remapping by shifted genome 
unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
PBMC_chrM_unmapped_bam=PBMC_possorted_chrM_and_unmapped.bam
nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_shifted.bam - &

# filtered unmapped reads 
samtools index PBMC_chrM_unmapped_bwa_shifted.bam
samtools view -@ 12 -b PBMC_chrM_unmapped_bwa_shifted.bam chrM > PBMC_chrM_unmapped_bwa_shifted2.bam
# only keep the DLoop region chrM:7000-10000 
samtools index PBMC_chrM_unmapped_bwa_shifted2.bam
samtools view -h -b PBMC_chrM_unmapped_bwa_shifted2.bam chrM:7000-10000 > PBMC_chrM_Dloop.bam

# Step3:  sort by CB 
all_bam=/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/PBMC_chrM_mapped.bam
Dloop_bam=/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/PBMC_chrM_Dloop.bam

samtools sort -@ 36 -t CB $all_bam -o PBMC_chrM_sorted_CB.bam
samtools sort -@ 36 -t CB $Dloop_bam -o PBMC_chrM_Dloop_sorted_CB.bam

# input a barcode,celltype annotation file 
cut -d, -f1 /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv | tr -d '"' > barcodes.txt

# Step5: split bam by barcode 
python split_bam.py -i PBMC_chrM_sorted_CB.bam -b barcodes.txt -o ./unshifted_bam/
python split_bam.py -i PBMC_chrM_Dloop_sorted_CB.bam -b barcodes.txt -o ./Dloop_bam/

# Step6: mutation calling for each barcode;

# nohup sh ./SNV_calling_v1.sh > SNV_Calling_outputv1.log 2>&1 &
# nohup sh ./SNV_calling_v2.sh > SNV_Calling_outputv2.log 2>&1 &

nohup sh ./SNV_calling.sh > SNV_Calling_output_allcell.log 2>&1 &


# 批量处理Step1-5：

# 1. sample_list.txt
# 格式：样本名 BAM文件路径 Barcode文件路径

PBMC_lib5_testppl /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv
PBMC_lib1 /data/R01/renwx5/MT_mution/outs/PBMC_lib_1/outs/atac_possorted_bam.bam /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation/Lib1_barcode_sample_annotation.csv
PBMC_lib3 /data/R01/renwx5/MT_mution/outs/PBMC_lib_3/outs/atac_possorted_bam.bam /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation/Lib3_barcode_sample_annotation.csv
PBMC_lib5 /data/R01/renwx5/MT_mution/outs/PBMC_lib_5/outs/atac_possorted_bam.bam /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation/Lib5_barcode_sample_annotation.csv

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
SHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
UNSHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"

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



nohup ./batch_process_mito.sh sample_list.txt > bam_processing.log 2>&1 &

# Step6: SNV calling (注意！！这一步不要多个样本一起跑！！！)

./process_mito_variants.sh <csv_file> <output_base> <unshifted_bam_base> <shifted_bam_base>


nohup sh ./SNV_calling.sh > SNV_Calling_output_allcell.log 2>&1 &


