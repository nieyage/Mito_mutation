#!/bin/bash

# 设置路径
csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/test.csv"
output_base="/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_v2"
unshifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/splitted_unshift"
shifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/splitted_shift"
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"

# 并行处理参数
PARALLEL_JOBS=10

# 质量过滤参数
MIN_MAPQ=30
MIN_BASEQ=30
MIN_VAR_FREQ=0.000001
MIN_READS2=3

# 合并mpileup的坐标参数
SHIFTED_REGION1_START=8570
SHIFTED_REGION1_END=9144
UNSHIFTED_REGION_START=576
UNSHIFTED_REGION_END=16024
SHIFTED_REGION2_START=8025
SHIFTED_REGION2_END=8569

# ==================== 初始化 ====================

# 创建输出目录
echo "创建输出目录: ${output_base}"
mkdir -p "${output_base}"
mkdir -p "${output_base}/logs"
echo "  unshifted BAM目录: ${unshifted_bam_base}"
echo "  shifted BAM目录: ${shifted_bam_base}"

# ==================== 函数定义 ====================

# 函数：从CSV文件中提取目标barcode（无临时文件版本）
get_target_barcodes() {
    awk -F, 'NR>1 && $1 != "" {
        gsub(/"/, "", $1);
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1);
        if ($1 != "") print $1
    }' "${csv_file}"
}

# 函数：高效内存处理单个barcode（零磁盘IO版本）
process_barcode_memory() {
    local barcode="$1"
    local output_prefix="${output_base}/${barcode}"
    local log_file="${output_base}/logs/${barcode}.log"
    
    echo "========== 开始处理细胞: ${barcode} =========="
    
    # 设置输入文件路径
    local unshifted_bam="${unshifted_bam_base}/${barcode}.bam"
    local shifted_bam="${shifted_bam_base}/${barcode}.bam"

    if [[ ! -f "${unshifted_bam}" ]]; then
        echo "错误: 找不到unshifted BAM文件: ${unshifted_bam}"
        return 1
    fi
    
    if [[ ! -f "${shifted_bam}" ]]; then
        echo "错误: 找不到shifted BAM文件: ${shifted_bam}"
        return 1
    fi
    
    # 步骤1: 使用命名管道处理，避免中间文件
    echo "步骤1-3: 使用管道直接处理BAM文件..."
    
    # 创建命名管道
    local unshifted_pipe shifted_pipe combined_pipe
    unshifted_pipe=$(mktemp -u /dev/shm/unshifted_XXXXX)
    shifted_pipe=$(mktemp -u /dev/shm/shifted_XXXXX)
    combined_pipe=$(mktemp -u /dev/shm/combined_XXXXX)
    
    mkfifo "${unshifted_pipe}" "${shifted_pipe}" "${combined_pipe}"
    
    # 1. 并行处理unshifted BAM（直接输出到管道）
    (
        # 使用内存排序，直接输出mpileup格式
        samtools view -u -q ${MIN_MAPQ} -F 4 "${unshifted_bam}" | \
        samtools sort -n -@ 2 -m 1G - | \
        samtools mpileup -Q ${MIN_BASEQ} -f "${unshifted_chrM_ref}" - 2>/dev/null | \
        awk -v start=${UNSHIFTED_REGION_START} -v end=${UNSHIFTED_REGION_END} \
            'BEGIN{OFS="\t"} $2>=start && $2<=end {print}' \
        > "${unshifted_pipe}"
    ) &
    
    # 2. 并行处理shifted BAM（直接输出到管道）
    (
        # 处理shifted BAM，分为两个区域
        samtools view -u -q ${MIN_MAPQ} -F 4 "${shifted_bam}" | \
        samtools sort -n -@ 2 -m 1G - | \
        samtools mpileup -Q ${MIN_BASEQ} -f "${shifted_chrM_ref}" - 2>/dev/null | \
        awk -v start1=${SHIFTED_REGION1_START} -v end1=${SHIFTED_REGION1_END} \
             -v start2=${SHIFTED_REGION2_START} -v end2=${SHIFTED_REGION2_END} \
            'BEGIN{OFS="\t"} 
             {
                 if ($2>=start1 && $2<=end1) {
                     $2 = $2 - start1 + 1;
                     print
                 } else if ($2>=start2 && $2<=end2) {
                     $2 = $2 - start2 + 15701;
                     print
                 }
             }' \
        > "${shifted_pipe}"
    ) &
    
    # 3. 合并两个管道的数据到最终管道
    (
        # 先读取shifted管道的数据（区域1）
        cat "${shifted_pipe}"
        # 再读取unshifted管道的数据（中间区域）
        cat "${unshifted_pipe}"
        # 再读取shifted管道的数据（区域2）- 需要特殊处理
    ) > "${combined_pipe}" &
    
    # 4. 从合并管道中读取并处理
    echo "步骤4: 变异检测..."
    
    # 创建临时内存文件用于varscan处理
    local varscan_temp
    varscan_temp=$(mktemp /dev/shm/varscan_XXXXX)
        
    cat "${combined_pipe}" > "${varscan_temp}"
    cat "${varscan_temp}" | varscan pileup2snp \
        --min-var-freq ${MIN_VAR_FREQ} \
        --min-reads2 ${MIN_READS2} - \
        > "${output_prefix}.snv" 2>/dev/null
    
    echo "步骤5: 提取pileup信息..."
    
    # 使用保存的临时文件提取pileup信息
    if [[ -s "${varscan_temp}" ]]; then
        "${pileup_script}" "${varscan_temp}" > "${output_prefix}.counts" 2>/dev/null
    fi
    
    # 清理所有命名管道和临时文件
    rm -f "${unshifted_pipe}" "${shifted_pipe}" "${combined_pipe}" #"${varscan_temp}"
    
    # 等待所有后台进程完成
    wait
    
    echo "========== 完成处理细胞: ${barcode} =========="
    return 0
}

# ==================== 主程序 ====================

# 从CSV文件中提取目标barcode
echo "从CSV文件中读取目标barcode..."
TARGET_BARCODES=($(get_target_barcodes "${csv_file}"))

echo "找到 ${#TARGET_BARCODES[@]} 个目标barcode"

if [[ ${#TARGET_BARCODES[@]} -eq 0 ]]; then
    echo "错误: 未找到有效的barcode"
    exit 1
fi

# 显示前几个barcode以验证
echo "前5个barcode示例:"
for i in {0..4}; do
    [[ $i -lt ${#TARGET_BARCODES[@]} ]] && echo "  ${TARGET_BARCODES[$i]}"
done

# 导出函数和变量
export -f process_barcode_memory
export output_base unshifted_bam_base shifted_bam_base
export unshifted_chrM_ref shifted_chrM_ref
export MIN_MAPQ MIN_BASEQ MIN_VAR_FREQ MIN_READS2
export SHIFTED_REGION1_START SHIFTED_REGION1_END
export UNSHIFTED_REGION_START UNSHIFTED_REGION_END
export SHIFTED_REGION2_START SHIFTED_REGION2_END
export pileup_script

# 使用GNU Parallel并行处理（优化参数）
echo "开始并行处理 ${#TARGET_BARCODES[@]} 个细胞..."
echo "并行任务数: ${PARALLEL_JOBS}"

# 设置进程限制和资源控制
ulimit -n 4096  # 增加文件描述符限制

# 使用高效的parallel参数
printf '%s\n' "${TARGET_BARCODES[@]}" | \
parallel -j ${PARALLEL_JOBS} \
    --progress \
    --joblog "${output_base}/parallel_joblog.txt" \
    --eta \
    --results "${output_base}/parallel_results" \
    "process_barcode_memory {} "

# 统计处理结果
echo "=================== 处理结果统计 ==================="

success_count=0
total_count=${#TARGET_BARCODES[@]}

for barcode in "${TARGET_BARCODES[@]}"; do
    exit_file="${output_base}/logs/${barcode}.exit"
    if [[ -f "${exit_file}" ]]; then
        exit_code=$(cat "${exit_file}")
        if [[ ${exit_code} -eq 0 ]]; then
            ((success_count++))
        fi
        rm -f "${exit_file}"
    fi
done

success_rate=$(echo "scale=2; ${success_count} * 100 / ${total_count}" | bc 2>/dev/null || echo "N/A")

echo "总细胞数: ${total_count}"
echo "成功处理: ${success_count}"
echo "失败处理: $((total_count - success_count))"
echo "成功率: ${success_rate}%"

# 生成最终报告
{
    echo "处理完成时间: $(date)"
    echo "CSV文件: ${csv_file}"
    echo "输出目录: ${output_base}"
    echo "总细胞数: ${total_count}"
    echo "成功处理: ${success_count}"
    echo "失败处理: $((total_count - success_count))"
    echo "并行任务数: ${PARALLEL_JOBS}"
    echo "处理参数: MAPQ=${MIN_MAPQ}, BASEQ=${MIN_BASEQ}"
    echo "变异检测: min-var-freq=${MIN_VAR_FREQ}, min-reads2=${MIN_READS2}"
} > "${output_base}/final_report.txt"

# 显示磁盘使用情况
echo "=================== 磁盘使用情况 ==================="
du -sh "${output_base}" 2>/dev/null || echo "无法统计磁盘使用"

# 显示输出文件示例
if [[ ${success_count} -gt 0 ]]; then
    echo "输出文件示例:"
    find "${output_base}" -name "*.snv" -type f -size +0 | head -3 | while read -r file; do
        base_name=$(basename "${file}" .snv)
        snv_count=$(tail -n +2 "${file}" 2>/dev/null | wc -l 2>/dev/null || echo 0)
        echo "  - ${base_name}.snv (${snv_count} 个变异)"
        echo "  - ${base_name}.counts"
    done
fi

echo "=================== 脚本执行完成 ==================="
echo "详细日志请查看: ${output_base}/logs/"
echo "最终报告: ${output_base}/final_report.txt"

exit 0
