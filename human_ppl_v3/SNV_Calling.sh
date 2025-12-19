#!/bin/bash

# 设置路径

csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"
output_base="/md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/masked_SNVcalling_percell"

# BAM文件基础路径
unshifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/splitted_unshift"
shifted_bam_base="/md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/splitted_shift"

# 参考基因组和工具路径

unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"
chrM_len="/md01/nieyg/ref/hard-mask/hg38_hard_masked/hg38_chrM.len"

# 并行处理参数
PARALLEL_JOBS=20

# 质量过滤参数
MIN_MAPQ=30
MIN_BASEQ=30
MIN_VAR_FREQ=0.000001
MIN_READS2=3

# 合并mpileup的坐标参数 
# 注意：这些参数可能需要根据您的实际数据调整
SHIFTED_REGION1_START=8570  # 偏移BAM中第一个区域起始位置
SHIFTED_REGION1_END=9144    # 偏移BAM中第一个区域结束位置
UNSHIFTED_REGION_START=576  # 原始BAM中间区域起始位置
UNSHIFTED_REGION_END=16024  # 原始BAM中间区域结束位置
SHIFTED_REGION2_START=8025  # 偏移BAM中第二个区域起始位置
SHIFTED_REGION2_END=8569    # 偏移BAM中第二个区域结束位置


# ==================== 初始化 ====================

# 创建输出目录
echo "创建输出目录: ${output_base}"
mkdir -p "${output_base}"
echo "  unshifted BAM目录: ${unshifted_bam_base}"
echo "  shifted BAM目录: ${shifted_bam_base}"

# ==================== 函数定义 ====================

# 函数：从CSV文件中提取目标barcode
get_target_barcodes() {
    local csv_file="$1"
    local target_barcodes=()
    
    # 创建临时文件
    local temp_file
    temp_file=$(mktemp)
    
    # 使用awk提取第一列的barcode，去除引号和首尾空格
    awk -F, 'NR>1 && $1 != "" {
        gsub(/"/, "", $1);
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1);
        if ($1 != "") {
            print $1
        }
    }' "${csv_file}" > "${temp_file}"
    
    # 从临时文件读取barcode
    while IFS= read -r barcode; do
        if [[ -n "${barcode}" ]]; then
            target_barcodes+=("${barcode}")
        fi
    done < "${temp_file}"
    
    # 清理临时文件
    rm -f "${temp_file}"
    printf '%s\n' "${target_barcodes[@]}"
}

# 函数：处理单个barcode
process_single_barcode() {
    local barcode="$1"
    
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
    
    # 设置输出文件前缀
    local output_prefix="${output_base}/${barcode}"
    
    # 步骤1: 处理unshifted BAM（原始坐标）
    echo "步骤1: 处理unshifted BAM..."
    
    # 1.1 排序
    samtools sort "${unshifted_bam}" -o "${output_prefix}_unshifted.sorted.bam"
    
    # 1.2 标记和去除重复reads
    java -Xmx4g -jar "${picard_tool}" \
        MarkDuplicates \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_unshifted.sorted.bam" \
        OUTPUT="${output_prefix}_unshifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_unshifted.metrics"
    
    # 1.3 生成mpileup文件
    samtools mpileup \
        -l "${chrM_len}" \
        -q ${MIN_MAPQ} \
        -Q ${MIN_BASEQ} \
        -f "${unshifted_chrM_ref}" \
        -x "${output_prefix}_unshifted.rmdup.bam" \
        > "${output_prefix}_unshifted.rmdup.mpileup"
    
    # 步骤2: 处理shifted BAM（偏移坐标）
    echo "步骤2: 处理shifted BAM..."
    
    # 2.1 排序
    samtools sort "${shifted_bam}" -o "${output_prefix}_shifted.sorted.bam"
    
    # 2.2 标记和去除重复reads
    java -Xmx4g -jar "${picard_tool}" \
        MarkDuplicates \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_shifted.sorted.bam" \
        OUTPUT="${output_prefix}_shifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_shifted.metrics"
    
    # 2.3 生成mpileup文件
    samtools mpileup \
        -l "${chrM_len}" \
        -q ${MIN_MAPQ} \
        -Q ${MIN_BASEQ} \
        -f "${shifted_chrM_ref}" \
        -x "${output_prefix}_shifted.rmdup.bam" | \
        sed 's/^chrM_shifted_8000_bp\t/chrM\t/' > "${output_prefix}_shifted.rmdup.mpileup"
    
    # 步骤3: 合并mpileup文件
    echo "步骤3: 合并mpileup文件..."
    
    local shifted_mpileup="${output_prefix}_shifted.rmdup.mpileup"
    local unshifted_mpileup="${output_prefix}_unshifted.rmdup.mpileup"
    local output_mpileup="${output_prefix}_combined.mpileup"
    
    # 3.1 从偏移mpileup中提取第一区域（对应原始坐标开始部分）
    awk -v start=${SHIFTED_REGION1_START} \
        -v end=${SHIFTED_REGION1_END} \
        -v new_start=1 \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             $2 = $2 - start + new_start;
             print
         }' "${shifted_mpileup}" > "${output_prefix}_temp1"
    
    # 3.2 从原始mpileup中提取中间区域
    awk -v start=${UNSHIFTED_REGION_START} \
        -v end=${UNSHIFTED_REGION_END} \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             print
         }' "${unshifted_mpileup}" > "${output_prefix}_temp2"
    
    # 3.3 从偏移mpileup中提取第二区域（对应原始坐标结束部分）
    awk -v start=${SHIFTED_REGION2_START} \
        -v end=${SHIFTED_REGION2_END} \
        -v new_start=15701 \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             $2 = $2 - start + new_start;
             print
         }' "${shifted_mpileup}" > "${output_prefix}_temp3"
    
    # 3.4 合并三个部分，构建完整的线粒体基因组mpileup
    cat "${output_prefix}_temp1" \
        "${output_prefix}_temp2" \
        "${output_prefix}_temp3" > "${output_mpileup}"
    
    # 清理临时文件
    rm -f "${output_prefix}_temp1" "${output_prefix}_temp2" "${output_prefix}_temp3"
    
    # 步骤4: 变异检测
    echo "步骤4: 变异检测..."
    varscan pileup2snp "${output_mpileup}" \
        --min-var-freq ${MIN_VAR_FREQ} \
        --min-reads2 ${MIN_READS2} \
        > "${output_prefix}.snv"
    
    # 步骤5: 提取pileup信息
    echo "步骤5: 提取pileup信息..."
    "${pileup_script}" "${output_mpileup}" > "${output_prefix}.counts"
    
    # 步骤6: 清理中间文件
    echo "步骤6: 清理中间文件..."
    rm -f \
        "${output_prefix}_unshifted.sorted.bam" \
        "${output_prefix}_shifted.sorted.bam" \
        "${output_prefix}_unshifted.rmdup.bam" \
        "${output_prefix}_shifted.rmdup.bam" \
        "${output_prefix}_unshifted.rmdup.mpileup" \
        "${output_prefix}_shifted.rmdup.mpileup"
    
    echo "========== 完成处理细胞: ${barcode} =========="
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
    if [[ $i -lt ${#TARGET_BARCODES[@]} ]]; then
        echo "  ${TARGET_BARCODES[$i]}"
    fi
done

# 导出函数和变量，以便在子shell中使用
export -f process_single_barcode
export output_base unshifted_bam_base shifted_bam_base
export unshifted_chrM_ref shifted_chrM_ref picard_tool pileup_script chrM_len
export MIN_MAPQ MIN_BASEQ MIN_VAR_FREQ MIN_READS2
export SHIFTED_REGION1_START SHIFTED_REGION1_END
export UNSHIFTED_REGION_START UNSHIFTED_REGION_END
export SHIFTED_REGION2_START SHIFTED_REGION2_END

# 使用GNU Parallel并行处理
echo "开始并行处理 ${#TARGET_BARCODES[@]} 个细胞..."
echo "并行任务数: ${PARALLEL_JOBS}"

printf '%s\n' "${TARGET_BARCODES[@]}" | \
parallel -j ${PARALLEL_JOBS} \
    --progress \
    --joblog "${output_base}/parallel_joblog.txt" \
    --eta \
    --results "${output_base}/parallel_results" \
    "process_single_barcode {}"

# 检查并行处理结果
if [[ $? -eq 0 ]]; then
    echo "所有barcode处理完成！"
    echo "结果保存在: ${output_base}"
    echo "日志文件: ${output_base}/parallel_joblog.txt"
    
    # 统计成功处理的细胞数量
    successful_cells=$(find "${output_base}" -name "*.snv" -type f | wc -l)
    echo "成功处理的细胞数量: ${successful_cells}"
    
    # 显示一些示例输出文件
    echo "示例输出文件:"
    find "${output_base}" -name "*.snv" -type f | head -3 | while read file; do
        echo "  - $(basename "${file}")"
    done
else
    echo "警告: 并行处理过程中出现错误"
    echo "请查看日志文件: ${output_base}/parallel_joblog.txt"
fi

# ==================== 结束 ====================
echo "脚本执行完成"
exit 0
