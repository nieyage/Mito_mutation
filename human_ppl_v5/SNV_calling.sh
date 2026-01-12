#!/bin/bash

csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"
output_base="/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_allcell"
unshifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/unshifted_bam"
shifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/Dloop_bam"
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"

# 并行处理参数
PARALLEL_JOBS=20

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
NEW_REGION1_START=1
NEW_REGION3_START=16025

# ==================== 初始化 ====================

# 创建输出目录
echo "创建输出目录: ${output_base}"
mkdir -p "${output_base}"

# ==================== 函数定义 ====================

# 函数：从CSV文件中提取目标barcode
get_target_barcodes() {
    local csv_file="$1"
    
    # 直接使用awk处理，避免创建临时文件
    awk -F, 'NR>1 && $1 != "" {
        gsub(/"/, "", $1);
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1);
        if ($1 != "") {
            print $1
        }
    }' "${csv_file}"
}

# 函数：在内存中处理mpileup区域提取
extract_mpileup_region() {
    local mpileup_input="$1"
    local start="$2"
    local end="$3"
    local new_start="$4"
    
    # 使用awk直接处理流数据，不写临时文件
    awk -v start="${start}" \
        -v end="${end}" \
        -v new_start="${new_start}" \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             $2 = $2 - start + new_start;
             print $0
         }' "${mpileup_input}"
}

process_single_barcode() {
    local barcode="$1"
    
    echo "========== 开始处理细胞: ${barcode} =========="
    
    local unshifted_bam="${unshifted_bam_base}/${barcode}.bam"
    local shifted_bam="${shifted_bam_base}/${barcode}.bam"
    local output_prefix="${output_base}/${barcode}"

    # 使用单个管道流处理所有步骤，避免中间文件
    echo "sort bam..."
    samtools sort "${unshifted_bam}" -o "${unshifted_bam_base}/${barcode}.sorted.bam"
    samtools sort "${shifted_bam}" -o "${shifted_bam_base}/${barcode}.sorted.bam"
    #rm "${unshifted_bam_base}/${barcode}.bam" "${shifted_bam_base}/${barcode}.bam"

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
export -f process_single_barcode extract_mpileup_region
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
        echo '处理完成: {} '
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

else
    echo "警告: 并行处理过程中出现错误"
    echo "请查看日志文件: ${output_base}/parallel_joblog.txt"
    echo "详细日志: ${output_base}/processing.log"
fi

# ==================== 结束 ====================
echo "脚本执行完成"
exit 0
