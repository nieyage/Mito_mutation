# Step1: make the ref 


#!/bin/bash

GENOME_FA="/md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa"  # 请修改为实际路径
OUTPUT_DIR="/md01/nieyg/ref/mito_ref/mm10"
SHIFT_BP=8000

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_command() {
    if ! command -v $1 &> /dev/null; then
        log_message "ERROR: $1 is not installed or not in PATH"
        exit 1
    fi
}

log_message "Starting chrM processing..."

check_command samtools
check_command awk

mkdir -p $OUTPUT_DIR

# 步骤1: 提取chrM序列
log_message "Step 1: Extracting chrM sequence..."
if [ -f "$GENOME_FA" ]; then
    samtools faidx $GENOME_FA chrM > $OUTPUT_DIR/mm10.chrM.fasta 2>/dev/null
    if [ $? -ne 0 ]; then
        log_message "Using alternative method to extract chrM..."
        awk 'BEGIN {RS=">"; ORS=""} /^>chrM[[:space:]]/ {print ">" $0}' $GENOME_FA > $OUTPUT_DIR/mm10.chrM.fasta
    fi
    samtools faidx $OUTPUT_DIR/mm10.chrM.fasta
else
    log_message "ERROR: Genome file not found: $GENOME_FA"
    exit 1
fi

ORIG_LENGTH=$(grep -v "^>" $OUTPUT_DIR/mm10.chrM.fasta | tr -d '\n' | wc -c)
log_message "Original chrM length: $ORIG_LENGTH bp"

# 步骤2: 位移8000bp
log_message "Step 2: Shifting chrM by $SHIFT_BP bp..."

# 创建临时工作目录
TEMP_DIR=$(mktemp -d)
cd $TEMP_DIR

# 将多行fasta转换为单行
grep -v "^>" $OUTPUT_DIR/mm10.chrM.fasta | tr -d '\n' > sequence.txt

# 调整位移大小（处理环状基因组）
ADJUSTED_SHIFT=$((SHIFT_BP % ORIG_LENGTH))
log_message "Adjusted shift (circular): $ADJUSTED_SHIFT bp"

# 执行位移
awk -v shift=$ADJUSTED_SHIFT -v len=$ORIG_LENGTH '
    {
        first_part = substr($0, shift + 1, len - shift);
        second_part = substr($0, 1, shift);
        printf("%s%s", first_part, second_part);
    }
' sequence.txt > shifted_sequence.txt

# 重建fasta格式
echo ">chrM_shifted_${SHIFT_BP}_bp" > header.txt
fold -w 80 shifted_sequence.txt > body.txt
cat header.txt body.txt > $OUTPUT_DIR/mm10.chrM.shifted_by_8000_bases.fasta

# 验证结果
SHIFTED_LENGTH=$(grep -v "^>" $OUTPUT_DIR/mm10.chrM.shifted_by_8000_bases.fasta | tr -d '\n' | wc -c)
if [ $ORIG_LENGTH -eq $SHIFTED_LENGTH ]; then
    log_message "✓ Shifted sequence length preserved: $SHIFTED_LENGTH bp"
else
    log_message "✗ ERROR: Length mismatch! Original: $ORIG_LENGTH, Shifted: $SHIFTED_LENGTH"
    exit 1
fi

# 创建shifted序列的索引
samtools faidx $OUTPUT_DIR/mm10.chrM.shifted_by_8000_bases.fasta

# 清理
cd /
rm -rf $TEMP_DIR

# ============================================
# 生成验证报告
# ============================================
log_message "Generating verification report..."
cat > $OUTPUT_DIR/chrM_processing_report.txt << EOF
chrM Sequence Processing Report
================================
Date: $(date)
Original genome: $GENOME_FA
Output directory: $OUTPUT_DIR

1. Original chrM sequence:
   File: mm10.chrM.fasta
   Length: $ORIG_LENGTH bp
   Index: mm10.chrM.fasta.fai (created)

2. Shifted chrM sequence:
   File: mm10.chrM.shifted_by_8000_bases.fasta
   Length: $SHIFTED_LENGTH bp
   Shift amount: $SHIFT_BP bp (adjusted to $ADJUSTED_SHIFT bp for circular genome)
   Index: mm10.chrM.shifted_by_8000_bases.fasta.fai (created)

3. Verification:
   Length preservation: $( [ $ORIG_LENGTH -eq $SHIFTED_LENGTH ] && echo "PASS" || echo "FAIL" )
   Files created successfully: $( [ -f "$OUTPUT_DIR/mm10.chrM.fasta" ] && [ -f "$OUTPUT_DIR/mm10.chrM.shifted_by_8000_bases.fasta" ] && echo "PASS" || echo "FAIL" )

EOF

# 显示前50个碱基对比
log_message "First 50 bases comparison:"
echo "Original:  $(grep -v '^>' $OUTPUT_DIR/mm10.chrM.fasta | tr -d '\n' | cut -c8001-8050)"
echo "Shifted:   $(grep -v '^>' $OUTPUT_DIR/mm10.chrM.shifted_by_8000_bases.fasta | tr -d '\n' | cut -c1-50)"

log_message "Processing complete!"
log_message "Output files:"
ls -lh $OUTPUT_DIR/mm10.chrM.*


# ref:https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
bwa index mm10.chrM.fasta
bwa index mm10.chrM.shifted_by_8000_bases.fasta
samtools faidx Homo_sapiens_assembly38.chrM.fasta
samtools faidx Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=mm10.chrM.fasta O=mm10.chrM.dict
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=mm10.chrM.shifted_by_8000_bases.fasta O=mm10.chrM.shifted_by_8000_bases.dict


# Step2: Get unmapped and chrM reads and Remapped with shifted chrM genome by bwa
input_bam=/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/BMMC_27m_ATAC/outs/atac_possorted_bam.bam
samtools view -b $input_bam chrM > BMMC_27m_chrM_mapped.bam
samtools view -b -f 4 $input_bam > BMMC_27m_unmapped.bam
samtools merge BMMC_27m_possorted_chrM_and_unmapped.bam BMMC_27m_unmapped.bam BMMC_27m_chrM_mapped.bam
samtools index BMMC_27m_possorted_chrM_and_unmapped.bam
rm BMMC_27m_unmapped.bam BMMC_27m_chrM_mapped.bam

# https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references

unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta
BMMC_27m_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_possorted_chrM_and_unmapped.bam
# remapping by shifted genome 
samtools collate -Oun128 $BMMC_27m_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $BMMC_27m_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o BMMC_27m_chrM_unmapped_bwa_shifted.bam - 

samtools index BMMC_27m_chrM_unmapped_bwa_shifted.bam

#(sleep 3600 && ./remap_chrM_shift.sh) &

nohup sh -c './remap_chrM_shift.sh' > output.log 2>&1 &

nohup samtools sort -@ 12 -t CB BMMC_27m_chrM_unmapped_bwa_shifted.bam -o BMMC_27m_chrM_unmapped_bwa_shifted_sorted_CB.bam & 
nohup samtools sort -@ 36 -t CB BMMC_27m_possorted_chrM_and_unmapped.bam -o BMMC_27m_chrM_unmapped_bwa_unshifted_sorted_CB.bam & 

# R get high quality cells 

tmp<- readRDS("/data/R04/liyh526/project/00_PLOG_aging/01_data/02_processed_data/BMMC_PBMC34_clear_celltype.rds")
data<- tmp@meta.data
data<-data[which(data$orig.ident=="27M_BMMC"),]
anno_data<- data.frame(barcode=data$atac_barcode,celltype=data$cell_type)
write.csv(anno_data,"/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/barcode_celltype.csv",row.names=F)
anno_data<- data.frame(barcode=data$gex_barcode,celltype=data$cell_type)
write.csv(anno_data,"/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv",row.names=F)

cut -d, -f1 /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv | tr -d '"' > barcodes.txt

# split bam by barcode 
nohup ./splitbam_unshifted.py -i BMMC_27m_chrM_unmapped_bwa_unshifted_sorted_CB.bam -b barcodes.txt > splitbam_unshifted_output.log 2>&1 &
nohup ./splitbam_shifted.py -i BMMC_27m_chrM_unmapped_bwa_shifted_sorted_CB.bam -b barcodes.txt > splitbam_shifted_output.log 2>&1 &

# for test


nohup sh ./SNV_Calling.sh > SNV_Calling_output3.log 2>&1 &
