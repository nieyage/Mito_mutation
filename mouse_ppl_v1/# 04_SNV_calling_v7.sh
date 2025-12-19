#!/bin/bash

############################################################
# 脚本名称：单细胞线粒体SNV检测并行处理脚本
# 描述：从单细胞测序数据中提取特定细胞的线粒体变异（SNV）
# 用途：检测线粒体异质性，进行单细胞线粒体突变分析
# 注意事项：需要GNU Parallel、samtools、picard、varscan等工具
############################################################

# ==================== 配置参数 ====================
# 注意：请根据实际情况修改以下路径

# 输入文件
CSV_FILE="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv"
# 输出目录
OUTPUT_BASE="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/SNV_calling_percell2"
 
# 工具路径
SUBSET_BAM_TOOL="/md01/nieyg/software/subset-bam_linux"
PICARD_TOOL="/md01/nieyg/software/picard-tools-2.4.1/picard.jar"

PILEUP_SCRIPT="/md01/jinxu/bin/pileup_inf_rj.pl"

# 参考基因组和长度文件
# 原始chrM参考基因组
UNSHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"
# 偏移8000bp的chrM参考基因组
SHIFTED_CHRM_REF="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta"
# chrM长度文件
CHRM_LEN="/md01/nieyg/ref/hard-mask/mm10_hard_masked/mm10_chrM.len"

# 输入BAM文件（包含所有细胞）
# 原始坐标的BAM文件
UNSHIFTED_BAM_ALL="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_possorted_chrM_and_unmapped.bam"
# 偏移坐标的BAM文件
SHIFTED_BAM_ALL="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_chrM_unmapped_bwa_shifted.bam"

# 并行处理参数
PARALLEL_JOBS=50  # 同时处理的细胞数量

# 质量过滤参数
MIN_MAPQ=30      # 最小比对质量
MIN_BASEQ=30     # 最小碱基质量
MIN_VAR_FREQ=0.000001  # 最小变异频率（0.0001%）
MIN_READS2=3     # 最小支持变异的reads数

# 合并mpileup的坐标参数（针对mm10线粒体基因组，总长16299bp）
SHIFTED_REGION1_START=8300  # 偏移BAM中第一个区域起始位置
SHIFTED_REGION1_END=8900    # 偏移BAM中第一个区域结束位置
UNSHIFTED_REGION_START=601  # 原始BAM中间区域起始位置
UNSHIFTED_REGION_END=15700  # 原始BAM中间区域结束位置
SHIFTED_REGION2_START=7701  # 偏移BAM中第二个区域起始位置
SHIFTED_REGION2_END=8299    # 偏移BAM中第二个区域结束位置

# ==================== 初始化 ====================

# 创建输出目录
echo "创建输出目录: ${OUTPUT_BASE}"
mkdir -p "${OUTPUT_BASE}"

# 检查必需工具
echo "检查必需工具..."
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "错误: 未找到工具 $1"
        echo "请安装或添加到PATH: $1"
        exit 1
    fi
}

check_tool "samtools"
check_tool "java"
check_tool "varscan"
check_tool "parallel"

# 检查输入文件
echo "检查输入文件..."
check_file() {
    if [[ ! -f "$1" ]]; then
        echo "错误: 文件不存在: $1"
        exit 1
    fi
}

check_file "${CSV_FILE}"
check_file "${UNSHIFTED_BAM_ALL}"
check_file "${SHIFTED_BAM_ALL}"
check_file "${UNSHIFTED_CHRM_REF}"
check_file "${SHIFTED_CHRM_REF}"
check_file "${CHRM_LEN}"
check_file "${PILEUP_SCRIPT}"

# ==================== 函数定义 ====================

# 函数：从CSV文件中提取目标barcode（处理带引号的情况）
get_target_barcodes() {
    local csv_file="$1"
    local target_barcodes=()
    
    # 注意：这里不要使用echo输出到stdout，否则会被捕获到数组中
    # 进度信息可以通过stderr输出，或者直接在主函数中输出
    
    # 创建临时文件来存储处理后的barcode
    local temp_file
    temp_file=$(mktemp)
    
    # 使用awk提取第一列的barcode，去除引号
    # 方法1：使用substr函数去除第一个和最后一个字符（如果被双引号包围）
    # 方法2：使用gsub函数替换所有双引号为空
    awk -F, 'NR>1 && $1 != "" {
        # 去除第一个字段的双引号
        gsub(/"/, "", $1);
        # 去除首尾空格
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
    
    # 设置输出文件前缀 - 确保barcode中没有特殊字符
    local output_prefix="${OUTPUT_BASE}/${barcode}"
    
    # 步骤1: 创建临时barcode文件
    echo "步骤1: 创建临时barcode文件..."
    local temp_barcode_file="${output_prefix}_barcode.txt"
    echo "${barcode}" > "${temp_barcode_file}"
    
    # 步骤2: 提取unshifted BAM（原始坐标）
    echo "步骤2: 提取unshifted BAM..."
    "${SUBSET_BAM_TOOL}" \
        --bam "${UNSHIFTED_BAM_ALL}" \
        --cell-barcodes "${temp_barcode_file}" \
        --cores 10 \
        --out-bam "${output_prefix}_unshifted.raw.bam"
    
    # 步骤3: 提取shifted BAM（偏移坐标）
    echo "步骤3: 提取shifted BAM..."
    "${SUBSET_BAM_TOOL}" \
        --bam "${SHIFTED_BAM_ALL}" \
        --cell-barcodes "${temp_barcode_file}" \
        --cores 10 \
        --out-bam "${output_prefix}_shifted.raw.bam"
    
    # 清理临时barcode文件
    rm -f "${temp_barcode_file}"
    
    # 检查提取的BAM文件是否为空
    if [[ ! -s "${output_prefix}_unshifted.raw.bam" ]] || [[ ! -s "${output_prefix}_shifted.raw.bam" ]]; then
        echo "警告: 细胞 ${barcode} 的BAM文件为空，跳过此细胞"
        rm -f "${output_prefix}_unshifted.raw.bam" "${output_prefix}_shifted.raw.bam"
        return 1
    fi
    
    # 步骤4: 处理unshifted BAM（原始坐标）
    echo "步骤4: 处理unshifted BAM..."
    
    # 4.1 排序
    samtools sort "${output_prefix}_unshifted.raw.bam" \
        -o "${output_prefix}_unshifted.sorted.bam"
    
    # 4.2 标记和去除重复reads
    java -Xmx4g -jar "${PICARD_TOOL}" \
        MarkDuplicates \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_unshifted.sorted.bam" \
        OUTPUT="${output_prefix}_unshifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_unshifted.metrics"
    
    # 4.3 生成mpileup文件
    samtools mpileup \
        -l "${CHRM_LEN}" \
        -q ${MIN_MAPQ} \
        -Q ${MIN_BASEQ} \
        -f "${UNSHIFTED_CHRM_REF}" \
        -x "${output_prefix}_unshifted.rmdup.bam" \
        > "${output_prefix}_unshifted.rmdup.mpileup"
    
    # 步骤5: 处理shifted BAM（偏移坐标）
    echo "步骤5: 处理shifted BAM..."
    
    # 5.1 排序
    samtools sort "${output_prefix}_shifted.raw.bam" \
        -o "${output_prefix}_shifted.sorted.bam"
    
    # 5.2 标记和去除重复reads
    java -Xmx4g -jar "${PICARD_TOOL}" \
        MarkDuplicates \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        INPUT="${output_prefix}_shifted.sorted.bam" \
        OUTPUT="${output_prefix}_shifted.rmdup.bam" \
        METRICS_FILE="${output_prefix}_shifted.metrics"
    
    # 5.3 生成mpileup文件
    samtools mpileup \
        -l "${CHRM_LEN}" \
        -q ${MIN_MAPQ} \
        -Q ${MIN_BASEQ} \
        -f "${SHIFTED_CHRM_REF}" \
        -x "${output_prefix}_shifted.rmdup.bam" \
        > "${output_prefix}_shifted.rmdup.mpileup"
    
    # 步骤6: 合并mpileup文件
    echo "步骤6: 合并mpileup文件..."
    
    local shifted_mpileup="${output_prefix}_shifted.rmdup.mpileup"
    local unshifted_mpileup="${output_prefix}_unshifted.rmdup.mpileup"
    local output_mpileup="${output_prefix}_combined.mpileup"
    
    # 6.1 从偏移mpileup中提取第一区域（对应原始坐标开始部分）
    awk -v start=${SHIFTED_REGION1_START} \
        -v end=${SHIFTED_REGION1_END} \
        -v new_start=1 \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             $2 = $2 - start + new_start;
             print
         }' "${shifted_mpileup}" > "${output_prefix}_temp1"
    
    # 6.2 从原始mpileup中提取中间区域
    awk -v start=${UNSHIFTED_REGION_START} \
        -v end=${UNSHIFTED_REGION_END} \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             print
         }' "${unshifted_mpileup}" > "${output_prefix}_temp2"
    
    # 6.3 从偏移mpileup中提取第二区域（对应原始坐标结束部分）
    awk -v start=${SHIFTED_REGION2_START} \
        -v end=${SHIFTED_REGION2_END} \
        -v new_start=15701 \
        'BEGIN{OFS="\t"} 
         $2>=start && $2<=end {
             $2 = $2 - start + new_start;
             print
         }' "${shifted_mpileup}" > "${output_prefix}_temp3"
    
    # 6.4 合并三个部分，构建完整的线粒体基因组mpileup
    cat "${output_prefix}_temp1" \
        "${output_prefix}_temp2" \
        "${output_prefix}_temp3" > "${output_mpileup}"
    
    # 清理临时文件
    rm -f "${output_prefix}_temp1" "${output_prefix}_temp2" "${output_prefix}_temp3"
    
    # 步骤7: 变异检测
    echo "步骤7: 变异检测..."
    varscan pileup2snp "${output_mpileup}" \
        --min-var-freq ${MIN_VAR_FREQ} \
        --min-reads2 ${MIN_READS2} \
        > "${output_prefix}.snv"
    
    # 步骤8: 提取pileup信息
    echo "步骤8: 提取pileup信息..."
    "${PILEUP_SCRIPT}" "${output_mpileup}" > "${output_prefix}.counts"
    
    # 步骤9: 清理中间文件
    echo "步骤9: 清理中间文件..."
    rm -f \
        "${output_prefix}_unshifted.raw.bam" \
        "${output_prefix}_shifted.raw.bam" \
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
echo "CSV文件: ${CSV_FILE}"
echo "提取barcode中..."

# 调用函数获取barcode数组
TARGET_BARCODES=($(get_target_barcodes "${CSV_FILE}"))
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
export OUTPUT_BASE UNSHIFTED_CHRM_REF SHIFTED_CHRM_REF PICARD_TOOL PILEUP_SCRIPT
export CHRM_LEN SUBSET_BAM_TOOL UNSHIFTED_BAM_ALL SHIFTED_BAM_ALL
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
    --joblog "${OUTPUT_BASE}/parallel_joblog.txt" \
    --eta \
    --results "${OUTPUT_BASE}/parallel_results" \
    "process_single_barcode {}"

# 检查并行处理结果
if [[ $? -eq 0 ]]; then
    echo "所有barcode处理完成！"
    echo "结果保存在: ${OUTPUT_BASE}"
    echo "日志文件: ${OUTPUT_BASE}/parallel_joblog.txt"
    
    # 统计成功处理的细胞数量
    successful_cells=$(find "${OUTPUT_BASE}" -name "*.snv" -type f | wc -l)
    echo "成功处理的细胞数量: ${successful_cells}"
    
    # 显示一些示例输出文件
    echo "示例输出文件:"
    find "${OUTPUT_BASE}" -name "*.snv" -type f | head -3 | while read file; do
        echo "  - $(basename "${file}")"
    done
else
    echo "警告: 并行处理过程中出现错误"
    echo "请查看日志文件: ${OUTPUT_BASE}/parallel_joblog.txt"
fi

# ==================== 结束 ====================
echo "脚本执行完成"
exit 0




nohup sh SNV_calling_v7_percell.sh > output.log 2>&1 &



