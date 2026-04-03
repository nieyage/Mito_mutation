#!/bin/bash

# ==================== 配置部分 ====================

# 基础路径配置
output_base="/md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/SNV_calling"
unshifted_bam_base="/md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/unshifted_bam"
shifted_bam_base="/md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/Dloop_bam"

# 工具与参考基因组路径
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
shifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
picard_tool="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup_script="/md01/nieyg/pipeline/mito_mutation/01_human/mpileup2readcounts.py"

# 并行处理参数 
PARALLEL_JOBS=7

# 过滤参数
MIN_MAPQ=30
MIN_BASEQ=30
MIN_VAR_FREQ=0.001
MIN_READS2=1

# 坐标转换参数
SHIFTED_REGION1_START=8570
SHIFTED_REGION1_END=9144
UNSHIFTED_REGION_START=576
UNSHIFTED_REGION_END=16024
SHIFTED_REGION2_START=8025
SHIFTED_REGION2_END=8569
NEW_REGION1_START=1
NEW_REGION3_START=16025

# ==================== 初始化 ====================

mkdir -p "${output_base}"
mkdir -p "${output_base}/logs"

# ==================== 函数定义 ====================

process_single_celltype() {
    local celltype_name="$1" # 传入的是类似 "CD14_Monocyte" 的名称
    
    echo "========== 开始处理细胞类型: ${celltype_name} =========="
    
    local unshifted_bam="${unshifted_bam_base}/${celltype_name}.bam"
    local shifted_bam="${shifted_bam_base}/${celltype_name}.bam"
    local output_prefix="${output_base}/${celltype_name}"

    # 检查输入文件是否存在
    if [[ ! -f "${unshifted_bam}" ]] || [[ ! -f "${shifted_bam}" ]]; then
        echo "错误: 找不到 ${celltype_name} 的 BAM 文件 (请检查 unshifted/shifted 路径)"
        return 1
    fi

    # 1. 排序 (如果你的 BAM 已经排过序，可以跳过此步骤)
    echo "[${celltype_name}] Sorting BAMs..."
    samtools sort -@ 4 "${unshifted_bam}" -o "${output_prefix}.unshifted.sorted.bam"
    samtools sort -@ 4 "${shifted_bam}" -o "${output_prefix}.shifted.sorted.bam"

    # 2. 流式合并 mpileup
    echo "[${celltype_name}] Generating combined mpileup..."
    {
        # 处理 shifted BAM - 区域 1
        java -Xmx8g -jar "${picard_tool}" MarkDuplicates \
            INPUT="${output_prefix}.shifted.sorted.bam" OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}.metrics" REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true 2>/dev/null | \
        samtools mpileup -q ${MIN_MAPQ} -Q ${MIN_BASEQ} -f "${shifted_chrM_ref}" -x - 2>/dev/null | \
        awk -v start=${SHIFTED_REGION1_START} -v end=${SHIFTED_REGION1_END} -v new_start=${NEW_REGION1_START} \
            'BEGIN{OFS="\t"} $2>=start && $2<=end {$2=$2-start+new_start; print $0}'
        
        # 处理 unshifted BAM - 中间区域
        java -Xmx8g -jar "${picard_tool}" MarkDuplicates \
            INPUT="${output_prefix}.unshifted.sorted.bam" OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}.metrics" REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true 2>/dev/null | \
        samtools mpileup -q ${MIN_MAPQ} -Q ${MIN_BASEQ} -f "${unshifted_chrM_ref}" -x - 2>/dev/null | \
        awk -v start=${UNSHIFTED_REGION_START} -v end=${UNSHIFTED_REGION_END} \
            'BEGIN{OFS="\t"} $2>=start && $2<=end {print $0}'

        # 处理 shifted BAM - 区域 2
        java -Xmx8g -jar "${picard_tool}" MarkDuplicates \
            INPUT="${output_prefix}.shifted.sorted.bam" OUTPUT=/dev/stdout \
            METRICS_FILE="${output_prefix}.metrics" REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true 2>/dev/null | \
        samtools mpileup -q ${MIN_MAPQ} -Q ${MIN_BASEQ} -f "${shifted_chrM_ref}" -x - 2>/dev/null | \
        awk -v start=${SHIFTED_REGION2_START} -v end=${SHIFTED_REGION2_END} -v new_start=${NEW_REGION3_START} \
            'BEGIN{OFS="\t"} $2>=start && $2<=end {$2=$2-start+new_start; print $0}'
    } > "${output_prefix}_combined.mpileup"

    # 3. 变异检测与计数提取
    echo "[${celltype_name}] Running VarScan and Count script..."
    varscan pileup2snp "${output_prefix}_combined.mpileup" \
        --min-var-freq ${MIN_VAR_FREQ} --min-reads2 ${MIN_READS2} > "${output_prefix}.snv" 2>/dev/null &
    
    python "${pileup_script}" "${output_prefix}_combined.mpileup" "${output_prefix}.counts" 2>/dev/null &
    
    wait

    # 4. 清理
    rm -f "${output_prefix}_combined.mpileup" "${output_prefix}.metrics" \
          "${output_prefix}.unshifted.sorted.bam" "${output_prefix}.shifted.sorted.bam"
    
    echo "========== 完成处理细胞类型: ${celltype_name} =========="
}

export -f process_single_celltype
export -f get_target_barcodes # 虽然不在此处用，但保留结构

# 将所有配置变量导出给 parallel
export output_base unshifted_bam_base shifted_bam_base unshifted_chrM_ref shifted_chrM_ref \
       picard_tool pileup_script MIN_MAPQ MIN_BASEQ MIN_VAR_FREQ MIN_READS2 \
       SHIFTED_REGION1_START SHIFTED_REGION1_END UNSHIFTED_REGION_START UNSHIFTED_REGION_END \
       SHIFTED_REGION2_START SHIFTED_REGION2_END NEW_REGION1_START NEW_REGION3_START

# ==================== 主程序 ====================

# 自动获取 celltype 列表：扫描 unshifted_bam_base 目录下的所有 .bam 文件名
echo "扫描输入目录获取细胞类型..."
CELLTYPES=($(ls "${unshifted_bam_base}"/*.bam | xargs -n 1 basename | sed 's/\.bam//'))

echo "待处理的细胞类型: ${CELLTYPES[@]}"
echo "总计: ${#CELLTYPES[@]} 个类型"

# 并行执行
printf '%s\n' "${CELLTYPES[@]}" | \
parallel -j ${PARALLEL_JOBS} \
    --progress \
    --joblog "${output_base}/parallel_celltype_joblog.txt" \
    "{
        echo 'Starting: {}'
        process_single_celltype {}
    }" 2>&1 | tee "${output_base}/celltype_processing.log"

echo "脚本执行完成。"
