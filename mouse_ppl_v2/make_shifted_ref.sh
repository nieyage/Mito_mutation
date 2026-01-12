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