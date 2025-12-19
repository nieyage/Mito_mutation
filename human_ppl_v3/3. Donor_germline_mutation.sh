# Step1: split bam by sample 
# Sample ID from Mitosort:/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv

#!/bin/bash

# 输入文件
csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"  # 请替换为实际的CSV文件路径
unshifted_bam="/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_unshifted.realign.bam"
shifted_bam="/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_shifted.realign.bam"
subset_bam_tool="/md01/nieyg/software/subset-bam_linux"

# 输出目录
output_dir="./donor_bams/"
mkdir -p "$output_dir"

# 临时目录存放barcode列表
barcode_dir="./barcode_lists/"
mkdir -p "$barcode_dir"

# 从CSV文件中提取每个Donor的barcode
echo "Extracting barcodes for each donor from CSV file..."

# 提取Donor1的barcode
awk -F',' '$2 == "\"Donor1\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor1_barcodes.txt"
echo "Donor1: $(wc -l < ${barcode_dir}donor1_barcodes.txt) barcodes"

# 提取Donor2的barcode
awk -F',' '$2 == "\"Donor2\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor2_barcodes.txt"
echo "Donor2: $(wc -l < ${barcode_dir}donor2_barcodes.txt) barcodes"

# 提取Donor3的barcode
awk -F',' '$2 == "\"Donor3\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor3_barcodes.txt"
echo "Donor3: $(wc -l < ${barcode_dir}donor3_barcodes.txt) barcodes"

# 提取Donor4的barcode
awk -F',' '$2 == "\"Donor4\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor4_barcodes.txt"
echo "Donor4: $(wc -l < ${barcode_dir}donor4_barcodes.txt) barcodes"

# 为每个Donor提取BAM文件
echo "Extracting BAM files for each donor..."

process_donor() {
    local donor=$1
    local barcode_file="${barcode_dir}${donor}_barcodes.txt"
    
    echo "Processing $donor..."
    
    # 提取unshifted BAM
    echo "  Extracting unshifted BAM for $donor..."
    "$subset_bam_tool" --bam "$unshifted_bam" --cell-barcodes "$barcode_file" --cores 1 --out-bam "${output_dir}${donor}_unshifted.bam"
    
    # 提取shifted BAM
    echo "  Extracting shifted BAM for $donor..."
    "$subset_bam_tool" --bam "$shifted_bam" --cell-barcodes "$barcode_file" --cores 1 --out-bam "${output_dir}${donor}_shifted.bam"
    
    echo "  Finished processing $donor"
}

# 处理所有Donor
process_donor "donor1"
process_donor "donor2"
process_donor "donor3"
process_donor "donor4"

echo "All Done!"
echo "Output files are in: $output_dir"


# Step2: Call germline mutation for each donor 

#!/bin/bash
# 设置变量
unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
chrM_len="/md01/nieyg/ref/hard-mask/hg38_hard_masked/hg38_chrM.len"
picard_jar="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup2snp_tool="varscan pileup2snp"  # 请确保这个工具在PATH中
pileup_inf_tool="/md01/jinxu/bin/pileup_inf_rj.pl"

# donor 列表
donors=("donor1" "donor2" "donor3" "donor4")

# 处理每个 donor
for donor in "${donors[@]}"; do
    echo "Processing $donor..."
    
    # 创建 donor 目录（如果不存在）
    mkdir -p "$donor"
    
    # 定义文件路径
    unshifted_bam="./donor_bams/${donor}_unshifted.bam"
    shifted_bam="./donor_bams/${donor}_shifted.bam"
    
    # 检查输入文件是否存在
    if [[ ! -f "$unshifted_bam" ]]; then
        echo "Error: Unshifted BAM file not found: $unshifted_bam"
        continue
    fi
    if [[ ! -f "$shifted_bam" ]]; then
        echo "Error: Shifted BAM file not found: $shifted_bam"
        continue
    fi
    
    # 1. 处理 unshifted BAM 文件
    echo "  Processing unshifted BAM..."
    
    # 排序
    samtools sort "$unshifted_bam" -o "${donor}/${donor}_unshifted.sorted.bam"
    
    # 标记重复
    java -Xmx4g -jar "$picard_jar" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${donor}/${donor}_unshifted.sorted.bam" \
        OUTPUT="${donor}/${donor}_unshifted.rmdup.bam" \
        METRICS_FILE="${donor}/${donor}_unshifted.metrics"
    
    # 生成 mpileup
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$unshifted_chrM_ref" \
        -x "${donor}/${donor}_unshifted.rmdup.bam" > "${donor}/${donor}_unshifted.rmdup.mpileup"
    
    # 2. 处理 shifted BAM 文件
    echo "  Processing shifted BAM..."
    
    # 排序
    samtools sort "$shifted_bam" -o "${donor}/${donor}_shifted.sorted.bam"
    
    # 标记重复
    java -Xmx4g -jar "$picard_jar" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${donor}/${donor}_shifted.sorted.bam" \
        OUTPUT="${donor}/${donor}_shifted.rmdup.bam" \
        METRICS_FILE="${donor}/${donor}_shifted.metrics"
    
    # 生成 mpileup
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$shifted_chrM_ref" \
        -x "${donor}/${donor}_shifted.rmdup.bam" > "${donor}/${donor}_shifted.rmdup.mpileup"
    
    # 3. 合并 mpileup 文件
    echo "  Merging mpileup files..."
    
    shifted_mpileup="${donor}/${donor}_shifted.rmdup.mpileup"
    unshifted_mpileup="${donor}/${donor}_unshifted.rmdup.mpileup"
    output_mpileup="${donor}/${donor}_combined.mpileup"
    
    # 检查输入文件是否存在
    if [[ ! -f "$shifted_mpileup" ]]; then
        echo "Error: shifted mpileup file not found: $shifted_mpileup"
        continue
    fi
    if [[ ! -f "$unshifted_mpileup" ]]; then
        echo "Error: unshifted mpileup file not found: $unshifted_mpileup"
        continue
    fi
    
    # 临时文件
    temp_part1="${donor}/temp_part1.mpileup"
    temp_part2="${donor}/temp_part2.mpileup"
    temp_part3="${donor}/temp_part3.mpileup"
    
    # 清理函数
    cleanup() {
        rm -f "$temp_part1" "$temp_part2" "$temp_part3"
    }
    trap cleanup EXIT
    
    # 第一部分：从shifted文件中提取8570-9144行，位置改为1-575
    awk -v start=8570 -v end=9144 -v new_start=1 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        new_pos = $2 - start + new_start;
        $2 = new_pos;
        print $0;
    }' "$shifted_mpileup" > "$temp_part1"
    
    # 第二部分：从unshifted文件中提取576-16024行，位置不变
    awk -v start=576 -v end=16024 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        print $0;
    }' "$unshifted_mpileup" > "$temp_part2"
    
    # 第三部分：从shifted文件中提取8025-8569行，位置改为16025-16569
    awk -v start=8025 -v end=8569 -v new_start=16025 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        new_pos = $2 - start + new_start;
        $2 = new_pos;
        print $0;
    }' "$shifted_mpileup" > "$temp_part3"
    
    # 合并三个部分
    cat "$temp_part1" "$temp_part2" "$temp_part3" > "$output_mpileup"
    
    # 清理临时文件
    # cleanup
    trap - EXIT
    
    # 验证合并结果
    total_lines=$(wc -l < "$output_mpileup")
    echo "  Combined mpileup lines: $total_lines"
    
    # 4. SNV calling 和计数
    echo "  Performing SNV calling and counting..."
    
    # SNV calling
    $pileup2snp_tool "$output_mpileup" --min-var-freq 0.01 --min-reads2 3 > "${donor}/${donor}.snv"
    
    # 计数
    $pileup_inf_tool "$output_mpileup" > "${donor}/${donor}.counts"
    
    echo "  Finished processing $donor"
    echo "----------------------------------------"
done

echo "All donors processed successfully!"


# for donors

library(dplyr)
donor1_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor1/donor1.snv", header = TRUE)
germline <- filter(donor1_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor1_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor2_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor2/donor2.snv", header = TRUE)
germline <- filter(donor2_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor2_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor3_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor3/donor3.snv", header = TRUE)
germline <- filter(donor3_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor3_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor4_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor4/donor4.snv", header = TRUE)
germline <- filter(donor4_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor4_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

