#!/bin/bash
# 脚本名称：process_mito_variants.sh
# 用途：并行处理线粒体SNV检测
# 用法：./process_mito_variants.sh <csv_file> <output_base> <unshifted_bam_base> <shifted_bam_base>

set -e  # 遇到错误时退出

# ==================== 参数检查 ====================
if [ $# -lt 4 ]; then
    echo "用法: $0 <csv_file> <output_base> <unshifted_bam_base> <shifted_bam_base>"
    echo "示例: $0 /path/to/human-mix-info.csv /path/to/output /path/to/unshifted_bam /path/to/Dloop_bam"
    echo ""
    echo "参数说明:"
    echo "  csv_file:          包含barcode信息的CSV文件"
    echo "  output_base:       输出结果的基础目录"
    echo "  unshifted_bam_base:未偏移BAM文件的基础目录"
    echo "  shifted_bam_base:  偏移BAM文件的基础目录"
    exit 1
fi

# 从命令行参数获取路径
csv_file="$1"
output_base="$2"
unshifted_bam_base="$3"
shifted_bam_base="$4"

# ==================== 其他配置参数 ====================
# 参考基因组文件
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta"

# 工具路径
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"

# 并行处理参数
PARALLEL_JOBS=10

# 质量过滤参数
MIN_MAPQ=30
MIN_BASEQ=30
MIN_VAR_FREQ=0.000001
MIN_READS2=3

# 合并mpileup的坐标参数 

SHIFTED_REGION1_START=8300  # 偏移BAM中第一个区域起始位置
SHIFTED_REGION1_END=8900    # 偏移BAM中第一个区域结束位置
UNSHIFTED_REGION_START=601  # 原始BAM中间区域起始位置
UNSHIFTED_REGION_END=15700  # 原始BAM中间区域结束位置
SHIFTED_REGION2_START=7701  # 偏移BAM中第二个区域起始位置
SHIFTED_REGION2_END=8299    # 偏移BAM中第二个区域结束位置
NEW_REGION1_START=1
NEW_REGION3_START=15701


# ==================== 参数验证 ====================
echo "=========================================="
echo "参数验证:"
echo "=========================================="
echo "CSV文件:           $csv_file"
echo "输出目录:          $output_base"
echo "未偏移BAM目录:     $unshifted_bam_base"
echo "偏移BAM目录:       $shifted_bam_base"
echo ""

# 检查必要文件是否存在
check_file_exists() {
    if [[ ! -f "$1" ]] && [[ ! -d "$1" ]]; then
        echo "错误: $2 '$1' 不存在或无法访问"
        exit 1
    fi
}

check_file_exists "$csv_file" "CSV文件"
check_file_exists "$unshifted_bam_base" "未偏移BAM目录"
check_file_exists "$shifted_bam_base" "偏移BAM目录"
check_file_exists "$unshifted_chrM_ref" "未偏移参考基因组"
check_file_exists "$shifted_chrM_ref" "偏移参考基因组"
check_file_exists "$picard_tool" "Picard工具"
check_file_exists "$pileup_script" "Pileup脚本"

# ==================== 初始化 ====================

# 创建输出目录
echo "创建输出目录: ${output_base}"
mkdir -p "${output_base}"

# ==================== 函数定义 ====================

# 函数：从CSV文件中提取目标barcode
get_target_barcodes() {
    local csv_file="$1"
    
    # 直接使用awk处理，避免创建临时文件
    awk -F, 'NR>=1 && $1 != "" {
        gsub(/"/, "", $1);
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1);
        if ($1 != "") {
            print $1
        }
    }' "${csv_file}"
}

process_single_barcode() {
    local barcode="$1"
    
    echo "========== 开始处理细胞: ${barcode} =========="
    
    local unshifted_bam="${unshifted_bam_base}/${barcode}.bam"
    local shifted_bam="${shifted_bam_base}/${barcode}.bam"
    local output_prefix="${output_base}/${barcode}"

    # 检查输入BAM文件是否存在
    if [[ ! -f "${unshifted_bam}" ]] || [[ ! -f "${shifted_bam}" ]]; then
        echo "警告: BAM文件不存在，跳过细胞 ${barcode}"
        echo "  unshifted: ${unshifted_bam}"
        echo "  shifted: ${shifted_bam}"
        return 1
    fi
    
    # 创建细胞特定的输出目录
    mkdir -p "$(dirname "${output_prefix}")"

    # 使用单个管道流处理所有步骤，避免中间文件
    echo "排序BAM文件..."
    samtools sort "${unshifted_bam}" -o "${unshifted_bam_base}/${barcode}.sorted.bam" 2>/dev/null
    samtools sort "${shifted_bam}" -o "${shifted_bam_base}/${barcode}.sorted.bam" 2>/dev/null
    
    # 检查排序后的文件
    if [[ ! -s "${unshifted_bam_base}/${barcode}.sorted.bam" ]] || [[ ! -s "${shifted_bam_base}/${barcode}.sorted.bam" ]]; then
        echo "错误: 排序后的BAM文件为空，跳过细胞 ${barcode}"
        return 1
    fi

    echo "步骤1-3: 流式处理并合并mpileup..."
    # 创建合并的mpileup文件，但使用流式处理减少内存压力
    {
        # 处理shifted BAM - 区域1
        java -Xmx4g -jar "${picard_tool}" \
            MarkDuplicates \
            CREATE_INDEX=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT \
            REMOVE_DUPLICATES=true \
            INPUT="${shifted_bam_base}/${barcode}.sorted.bam" \
            OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}_unshifted.metrics" 2>/dev/null \
        | samtools view -b - 2>/dev/null \
        | samtools mpileup \
            -q ${MIN_MAPQ} \
            -Q ${MIN_BASEQ} \
            -f "${shifted_chrM_ref}" \
            -x - 2>/dev/null | \
        awk -v start=${SHIFTED_REGION1_START} \
            -v end=${SHIFTED_REGION1_END} \
            -v new_start=${NEW_REGION1_START} \
            'BEGIN{OFS="\t"} 
             $2>=start && $2<=end {
                 $2 = $2-start+new_start
                 print $0
             }'
        
        # 处理unshifted BAM - 中间区域
        java -Xmx4g -jar "${picard_tool}" \
            MarkDuplicates \
            CREATE_INDEX=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT \
            REMOVE_DUPLICATES=true \
            INPUT="${unshifted_bam_base}/${barcode}.sorted.bam" \
            OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}_unshifted.metrics" 2>/dev/null \
        | samtools view -b - 2>/dev/null \
        | samtools mpileup \
            -q ${MIN_MAPQ} \
            -Q ${MIN_BASEQ} \
            -f "${unshifted_chrM_ref}" \
            -x - 2>/dev/null | \
        awk -v start=${UNSHIFTED_REGION_START} \
            -v end=${UNSHIFTED_REGION_END} \
            'BEGIN{OFS="\t"} 
             $2>=start && $2<=end {
                 print $0
             }'

        
        # 处理shifted BAM - 区域2
        java -Xmx4g -jar "${picard_tool}" \
            MarkDuplicates \
            CREATE_INDEX=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT \
            REMOVE_DUPLICATES=true \
            INPUT="${shifted_bam_base}/${barcode}.sorted.bam" \
            OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}_unshifted.metrics" 2>/dev/null \
        | samtools view -b - 2>/dev/null \
        | samtools mpileup \
            -q ${MIN_MAPQ} \
            -Q ${MIN_BASEQ} \
            -f "${shifted_chrM_ref}" \
            -x - 2>/dev/null | \
        awk -v start=${SHIFTED_REGION2_START} \
            -v end=${SHIFTED_REGION2_END} \
            -v new_start=${NEW_REGION3_START} \
            'BEGIN{OFS="\t"} 
             $2>=start && $2<=end {
                 $2 = $2 - start + new_start
                 print $0
             }'
    } > "${output_prefix}_combined.mpileup"
    
    # 检查文件是否生成
    if [[ ! -s "${output_prefix}_combined.mpileup" ]]; then
        echo "错误: 合并的mpileup文件为空" >&2
        return 1
    fi
    
    echo "步骤4-5: 并行处理最终输出..."
    
    # 并行处理最终输出
    (
        # 变异检测
        varscan pileup2snp "${output_prefix}_combined.mpileup" \
            --min-var-freq ${MIN_VAR_FREQ} \
            --min-reads2 ${MIN_READS2} \
            > "${output_prefix}.snv" 2>/dev/null
    ) &
    
    (
        # 提取pileup信息
        "${pileup_script}" "${output_prefix}_combined.mpileup" > "${output_prefix}.counts" 2>/dev/null
    ) &
    
    wait
    
    echo "步骤6: 清理中间文件..."
    # 删除中间文件
    rm -f "${output_prefix}_combined.mpileup" "${output_prefix}_unshifted.metrics"
    rm -f "${unshifted_bam_base}/${barcode}.sorted.bam" "${shifted_bam_base}/${barcode}.sorted.bam"
    
    echo "========== 完成处理细胞: ${barcode} =========="
    return 0
}

# ==================== 主程序 ====================

# 从CSV文件中提取目标barcode
echo "从CSV文件中读取目标barcode..."
TARGET_BARCODES=($(get_target_barcodes "${csv_file}"))
TARGET_COUNT=${#TARGET_BARCODES[@]}
echo "找到 ${TARGET_COUNT} 个目标barcode"

if [[ ${TARGET_COUNT} -eq 0 ]]; then
    echo "错误: 未找到有效的barcode"
    exit 1
fi

# 显示前几个barcode以验证
echo "前5个barcode示例:"
for ((i=0; i<5 && i<TARGET_COUNT; i++)); do
    echo "  ${TARGET_BARCODES[$i]}"
done

PROCESS_FUNCTION="process_single_barcode"

# 导出函数和变量
export -f process_single_barcode
export output_base unshifted_bam_base shifted_bam_base
export unshifted_chrM_ref shifted_chrM_ref picard_tool pileup_script
export MIN_MAPQ MIN_BASEQ MIN_VAR_FREQ MIN_READS2
export SHIFTED_REGION1_START SHIFTED_REGION1_END
export UNSHIFTED_REGION_START UNSHIFTED_REGION_END
export NEW_REGION1_START NEW_REGION3_START
export SHIFTED_REGION2_START SHIFTED_REGION2_END
export PROCESS_FUNCTION

# 使用GNU Parallel并行处理
echo "开始并行处理 ${TARGET_COUNT} 个细胞..."
echo "并行任务数: ${PARALLEL_JOBS}"

# 创建日志目录
mkdir -p "${output_base}/logs"

# 使用更高效的并行处理
printf '%s\n' "${TARGET_BARCODES[@]}" | \
parallel -j ${PARALLEL_JOBS} \
    --progress \
    --joblog "${output_base}/parallel_joblog.txt" \
    --results "${output_base}/parallel_results" \
    --eta \
    --halt soon,fail=10% \
    "{
        echo '处理开始: {}'
        start_time=\$(date +%s)
        
        # 执行处理函数
        if ${PROCESS_FUNCTION} {} 2>&1; then
            status='成功'
        else
            status='失败'
        fi
        
        end_time=\$(date +%s)
        duration=\$((end_time - start_time))
        echo '处理完成: {} (耗时: \${duration}秒)'
    }" 2>&1 | tee "${output_base}/processing.log"

# 检查并行处理结果
if [[ $? -eq 0 ]]; then
    echo "所有barcode处理完成！"
    echo "结果保存在: ${output_base}"
    echo "日志文件: ${output_base}/parallel_joblog.txt"
    
    # 统计成功处理的细胞数量
    successful_cells=$(find "${output_base}" -name "*.snv" -type f 2>/dev/null | wc -l)
    failed_cells=$((TARGET_COUNT - successful_cells))
    
    echo "成功处理的细胞数量: ${successful_cells}"
    echo "失败的细胞数量: ${failed_cells}"
    
    # 生成摘要报告
    echo "生成处理摘要..."
    cat > "${output_base}/summary.txt" << EOF
处理摘要
==========
处理时间: $(date)
CSV文件: ${csv_file}
输出目录: ${output_base}
未偏移BAM目录: ${unshifted_bam_base}
偏移BAM目录: ${shifted_bam_base}

统计结果:
总barcode数: ${TARGET_COUNT}
成功处理: ${successful_cells}
失败处理: ${failed_cells}
成功率: $(awk "BEGIN {printf \"%.1f%%\", ${successful_cells}/${TARGET_COUNT}*100}")

输出文件:
- 每个细胞生成两个文件: <barcode>.snv 和 <barcode>.counts
- 并行处理日志: parallel_joblog.txt
- 详细处理日志: processing.log
- 各细胞详细结果: parallel_results/

注意: 如果失败细胞较多，请检查日志文件了解具体原因。
EOF
    
    echo "摘要报告已保存到: ${output_base}/summary.txt"

else
    echo "警告: 并行处理过程中出现错误"
    echo "请查看日志文件: ${output_base}/parallel_joblog.txt"
    echo "详细日志: ${output_base}/processing.log"
fi

# ==================== 结束 ====================
echo "脚本执行完成"
exit 0