#!/bin/bash

# 设置路径
csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"
output_base="/md01/nieyg/project/mito_mutation/01_pipeline/03_split_bam/masked_SNVcalling_percell"
subset_bam_tool="/md01/nieyg/software/subset-bam_linux"
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"
chrM_len="/md01/nieyg/ref/hard-mask/hg38_hard_masked/hg38_chrM.len"

# 原始大BAM文件路径
unshifted_bam_all="PBMC_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_possorted_chrM_and_unmapped.bam"
shifted_bam_all="/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_shifted.realign.bam"

# 创建输出目录
mkdir -p "$output_base"

# 从CSV文件中提取目标barcode
get_target_barcodes() {
    local csv_file="$1"
    local target_barcodes=()
    
    # 使用awk提取第一列的barcode
    while IFS= read -r barcode; do
        if [[ -n "$barcode" ]]; then
            target_barcodes+=("$barcode")
        fi
    done < <(awk -F, 'NR>0 && $1 != "" {print $1}' "$csv_file")
    
    printf '%s\n' "${target_barcodes[@]}"
}

echo "从CSV文件中读取目标barcode..."
target_barcodes=($(get_target_barcodes "$csv_file"))
echo "找到 ${#target_barcodes[@]} 个目标barcode"

# 创建处理单个barcode的函数
process_barcode() {
    local barcode="$1"
    local output_base="$2"
    local unshifted_chrM_ref="$3"
    local shifted_chrM_ref="$4"
    local picard_tool="$5"
    local pileup_script="$6"
    local chrM_len="$7"
    local subset_bam_tool="$8"
    local unshifted_bam_all="$9"
    local shifted_bam_all="${10}"
    
    # 设置输出文件前缀
    local output_prefix="$output_base/${barcode}"
    
    echo "开始处理 $barcode"
    
    # 创建临时barcode文件
    local temp_barcode_file="${output_prefix}_barcode.txt"
    echo "$barcode" > "$temp_barcode_file"
    
    # 使用subset-bam工具提取unshifted BAM
    echo "提取 $barcode 的unshifted BAM..."
    "$subset_bam_tool" --bam "$unshifted_bam_all" \
        --cell-barcodes "$temp_barcode_file" \
        --cores 1 \
        --out-bam "${output_prefix}_unshifted.raw.bam"
    
    # 使用subset-bam工具提取shifted BAM
    echo "提取 $barcode 的shifted BAM..."
    "$subset_bam_tool" --bam "$shifted_bam_all" \
        --cell-barcodes "$temp_barcode_file" \
        --cores 1 \
        --out-bam "${output_prefix}_shifted.raw.bam"
    
    # 清理临时barcode文件
    rm -f "$temp_barcode_file"
    
    # 检查提取的BAM文件是否为空
    if [[ ! -s "${output_prefix}_unshifted.raw.bam" ]] || [[ ! -s "${output_prefix}_shifted.raw.bam" ]]; then
        echo "警告: $barcode 的BAM文件为空，跳过"
        rm -f "${output_prefix}_unshifted.raw.bam" "${output_prefix}_shifted.raw.bam"
        return 1
    fi
    
    # 处理unshifted BAM
    samtools sort "${output_prefix}_unshifted.raw.bam" -o "${output_prefix}_unshifted.sorted.bam"
    
    java -Xmx4g -jar "$picard_tool" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_unshifted.sorted.bam" \
        OUTPUT="${output_prefix}_unshifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_unshifted.metrics"
    
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$unshifted_chrM_ref" \
        -x "${output_prefix}_unshifted.rmdup.bam" > "${output_prefix}_unshifted.rmdup.mpileup"
    
    # 处理shifted BAM
    samtools sort "${output_prefix}_shifted.raw.bam" -o "${output_prefix}_shifted.sorted.bam"
    
    java -Xmx4g -jar "$picard_tool" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_shifted.sorted.bam" \
        OUTPUT="${output_prefix}_shifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_shifted.metrics"
    
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$shifted_chrM_ref" \
        -x "${output_prefix}_shifted.rmdup.bam" > "${output_prefix}_shifted.rmdup.mpileup"
    
    # 合并mpileup文件
    local shifted_mpileup="${output_prefix}_shifted.rmdup.mpileup"
    local unshifted_mpileup="${output_prefix}_unshifted.rmdup.mpileup"
    local output_mpileup="${output_prefix}_combined.mpileup"
    
    awk -v start=8570 -v end=9144 -v new_start=1 'BEGIN{OFS="\t"} $2>=start&&$2<=end{$2=$2-start+new_start;print}' "$shifted_mpileup" > "${output_prefix}_temp1"
    awk -v start=576 -v end=16024 'BEGIN{OFS="\t"} $2>=start&&$2<=end{print}' "$unshifted_mpileup" > "${output_prefix}_temp2"
    awk -v start=8025 -v end=8569 -v new_start=16025 'BEGIN{OFS="\t"} $2>=start&&$2<=end{$2=$2-start+new_start;print}' "$shifted_mpileup" > "${output_prefix}_temp3"
    cat "${output_prefix}_temp1" "${output_prefix}_temp2" "${output_prefix}_temp3" > "$output_mpileup"
    rm -f "${output_prefix}_temp1" "${output_prefix}_temp2" "${output_prefix}_temp3"
    
    # 运行变异检测
    varscan pileup2snp "$output_mpileup" --min-var-freq 0.000001 --min-reads2 3 > "${output_prefix}.snv"
    
    # 运行pileup信息提取脚本
    "$pileup_script" "$output_mpileup" > "${output_prefix}.counts"
    
    # 清理中间文件以节省空间
    rm -f "${output_prefix}_unshifted.raw.bam" "${output_prefix}_shifted.raw.bam" \
          "${output_prefix}_unshifted.sorted.bam" "${output_prefix}_shifted.sorted.bam" \
          "${output_prefix}_unshifted.rmdup.bam" "${output_prefix}_shifted.rmdup.bam" \
          "${output_prefix}_unshifted.rmdup.mpileup" "${output_prefix}_shifted.rmdup.mpileup"
    
    echo "完成处理 $barcode"
}

# 导出函数，以便在子shell中使用
export -f process_barcode
export output_base unshifted_chrM_ref shifted_chrM_ref picard_tool pileup_script chrM_len subset_bam_tool unshifted_bam_all shifted_bam_all

# 使用GNU Parallel并行处理，同时运行50个任务
printf '%s\n' "${target_barcodes[@]}" | \
parallel -j 50 --progress --joblog "$output_base/parallel_joblog.txt" \
"process_barcode {} $output_base $unshifted_chrM_ref $shifted_chrM_ref $picard_tool $pileup_script $chrM_len $subset_bam_tool $unshifted_bam_all $shifted_bam_all"

echo "所有barcode处理完成！结果保存在: $output_base"

