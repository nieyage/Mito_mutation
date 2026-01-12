# Call germline mutation 

#!/bin/bash

# 设置路径
csv_file="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/barcodes.txt"
output_base="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/04_germline"
# BAM文件基础路径
unshifted_bam="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_possorted_chrM_and_unmapped.bam"
shifted_bam="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_chrM_unmapped_bwa_shifted.bam"

# 参考基因组和工具路径
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta"
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/jinxu/bin/pileup_inf_rj.pl"
chrM_len="/md01/nieyg/ref/hard-mask/mm10_hard_masked/mm10_chrM.len"

# 并行处理参数
PARALLEL_JOBS=20
# 质量过滤参数
MIN_MAPQ=30
MIN_BASEQ=30
MIN_VAR_FREQ=0.000001
MIN_READS2=3

# 合并mpileup的坐标参数（针对mm10线粒体基因组，总长16299bp）
# 注意：这些参数可能需要根据您的实际数据调整
SHIFTED_REGION1_START=8300  # 偏移BAM中第一个区域起始位置
SHIFTED_REGION1_END=8900    # 偏移BAM中第一个区域结束位置
UNSHIFTED_REGION_START=601  # 原始BAM中间区域起始位置
UNSHIFTED_REGION_END=15700  # 原始BAM中间区域结束位置
SHIFTED_REGION2_START=7701  # 偏移BAM中第二个区域起始位置
SHIFTED_REGION2_END=8299    # 偏移BAM中第二个区域结束位置
echo "创建输出目录: ${output_base}"
mkdir -p "${output_base}"
echo "  unshifted BAM目录: ${unshifted_bam}"
echo "  shifted BAM目录: ${shifted_bam}"

    # 设置输出文件前缀
    local output_prefix="BMMC_27m_germline"

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

    echo "========== 完成处理germline mutation calling =========="





