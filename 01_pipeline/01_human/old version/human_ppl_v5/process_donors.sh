#!/bin/bash
# 文件名: process_donors.sh

# 默认参数
BASE_DIR="/md01/nieyg/project/mito_mutation/01_pipeline/10_v5"
SCRIPTS_DIR="/md01/nieyg/pipeline/mito_mutation/01_human"
OUTPUT_DIR="${BASE_DIR}/donor_analysis"
THREADS=16
BATCH_SIZE=1000
WORKERS=10

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        -s|--scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        -w|--workers)
            WORKERS="$2"
            shift 2
            ;;
        -h|--help)
            echo "用法: $0 [选项]"
            echo "选项:"
            echo "  -b, --base-dir     基础目录 (默认: $BASE_DIR)"
            echo "  -s, --scripts-dir  脚本目录 (默认: $SCRIPTS_DIR)"
            echo "  -o, --output-dir   输出目录 (默认: $OUTPUT_DIR)"
            echo "  -t, --threads      线程数 (默认: $THREADS)"
            echo "      --batch-size   批次大小 (默认: $BATCH_SIZE)"
            echo "  -w, --workers      工作进程数 (默认: $WORKERS)"
            echo "  -h, --help         显示帮助信息"
            exit 0
            ;;
        *)
            echo "未知选项: $1"
            exit 1
            ;;
    esac
done

INPUT_DIR="${BASE_DIR}/"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 获取所有Donor目录
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "脚本目录: $SCRIPTS_DIR"
echo -e "\n正在扫描Donor目录..."

# 使用数组存储Donor目录
mapfile -t DONOR_DIRS < <(find "$INPUT_DIR" -maxdepth 1 -type d -name "Donor*" | sort)

if [ ${#DONOR_DIRS[@]} -eq 0 ]; then
    echo "错误：在 $INPUT_DIR 中未找到Donor目录"
    echo "尝试查找其他模式..."
    mapfile -t DONOR_DIRS < <(find "$INPUT_DIR" -maxdepth 1 -type d | grep -v "^$INPUT_DIR$" | sort)
    
    if [ ${#DONOR_DIRS[@]} -eq 0 ]; then
        echo "错误：未找到任何子目录"
        exit 1
    fi
fi

echo "找到 ${#DONOR_DIRS[@]} 个目录:"
for dir in "${DONOR_DIRS[@]}"; do
    donor_name=$(basename "$dir")
    file_count=$(find "$dir" -name "*.counts" -o -name "*.snv" 2>/dev/null | wc -l)
    echo "  - $donor_name ($file_count 个文件)"
done

echo -e "\n开始处理..."

# 处理日志
LOG_FILE="${OUTPUT_DIR}/processing_log_$(date +%Y%m%d_%H%M%S).txt"
{
echo "处理开始时间: $(date)"
echo "参数设置:"
echo "  基础目录: $BASE_DIR"
echo "  脚本目录: $SCRIPTS_DIR"
echo "  输出目录: $OUTPUT_DIR"
echo "  线程数: $THREADS"
echo "  批次大小: $BATCH_SIZE"
echo "  工作进程数: $WORKERS"
echo "========================================"
} | tee "$LOG_FILE"

# 循环处理每个Donor目录
for donor_dir in "${DONOR_DIRS[@]}"; do
    donor_name=$(basename "$donor_dir")
    donor_output="${OUTPUT_DIR}/${donor_name}"
    
    {
    echo -e "\n========================================"
    echo "开始处理: $donor_name"
    echo "输入目录: $donor_dir"
    echo "输出目录: $donor_output"
    echo "开始时间: $(date)"
    echo "========================================"
    } | tee -a "$LOG_FILE"
    
    # 创建Donor输出目录
    mkdir -p "$donor_output"
    
    # 切换到Donor输出目录
    cd "$donor_output" || {
        echo "无法切换到目录: $donor_output" | tee -a "$LOG_FILE"
        continue
    }
    
    # 步骤1: 提取coverage
    {
    echo "步骤1: 运行 extract_coverage.py"
    echo "命令: python ${SCRIPTS_DIR}/extract_coverage.py -i \"$donor_dir\" -o \"${donor_name}_barcode_coverage.tsv.gz\" -p $THREADS -b $BATCH_SIZE"
    } | tee -a "$LOG_FILE"
    
    python "${SCRIPTS_DIR}/extract_coverage.py" \
        -i "$donor_dir" \
        -o "${donor_name}_barcode_coverage.tsv.gz" \
        -p "$THREADS" \
        -b "$BATCH_SIZE" 2>&1 | tee -a "$LOG_FILE"
    
    STEP1_EXIT=$?
    
    # 步骤2: SNV汇总
    {
    echo -e "\n步骤2: 运行 SNV_summary_v2.py"
    echo "命令: python ${SCRIPTS_DIR}/SNV_summary_v2.py -i \"$donor_dir\" -o \"${donor_name}_snv.tsv\" -v"
    } | tee -a "$LOG_FILE"
    
    python "${SCRIPTS_DIR}/SNV_summary_v2.py" \
        -i "$donor_dir" \
        -o "${donor_name}_snv.tsv" \
        -v 2>&1 | tee -a "$LOG_FILE"
    
    STEP2_EXIT=$?
    
    # 压缩snv文件
    if [ $STEP2_EXIT -eq 0 ] && [ -f "${donor_name}_snv.tsv" ]; then
        echo "压缩snv文件..." | tee -a "$LOG_FILE"
        gzip -f "${donor_name}_snv.tsv"
        snv_file="${donor_name}_snv.tsv.gz"
    elif [ -f "${donor_name}_snv.tsv.gz" ]; then
        snv_file="${donor_name}_snv.tsv.gz"
        echo "使用已压缩的snv文件: $snv_file" | tee -a "$LOG_FILE"
    else
        echo "错误: 找不到snv文件，跳过后续步骤" | tee -a "$LOG_FILE"
        continue
    fi
    
    # 步骤3: 过滤突变
    {
    echo -e "\n步骤3: 运行 filter_mutations.py"
    echo "命令: python3 ${SCRIPTS_DIR}/filter_mutations.py -i \"$snv_file\" -o ./ --germline-output \"${donor_name}_germline.tsv\" --somatic-output \"${donor_name}_somatic.tsv\" --strand-min 0.3 --strand-max 0.7 --min-depth 10 --alt-ratio-threshold 0.90 --no-compress"
    } | tee -a "$LOG_FILE"
    
    python3 "${SCRIPTS_DIR}/filter_mutations.py" \
        -i "$snv_file" \
        -o ./ \
        --germline-output "${donor_name}_germline.tsv" \
        --somatic-output "${donor_name}_somatic.tsv" \
        --strand-min 0.3 \
        --strand-max 0.7 \
        --min-depth 10 \
        --alt-ratio-threshold 0.90 \
        --no-compress 2>&1 | tee -a "$LOG_FILE"
    
    STEP3_EXIT=$?
    
    # 确定somatic文件
    somatic_file=""
    for possible_file in \
        "${donor_name}_somatic_withoutstrand.tsv" \
        "somatic_withoutstrand.tsv" \
        "${donor_name}_somatic.tsv" \
        "somatic.tsv"; do
        if [ -f "$possible_file" ]; then
            somatic_file="$possible_file"
            break
        fi
    done
    
    if [ -z "$somatic_file" ]; then
        echo "警告: 未找到somatic文件，跳过创建变异矩阵" | tee -a "$LOG_FILE"
    else
        # 步骤4: 创建变异矩阵
        {
        echo -e "\n步骤4: 运行 create_variant_matrix_v2.py"
        echo "命令: python3 ${SCRIPTS_DIR}/create_variant_matrix_v2.py -s \"$somatic_file\" -d \"$donor_dir\" -o \"${donor_name}_variant_sparse_matrix_withoutstrand.tsv.gz\" --workers $WORKERS"
        } | tee -a "$LOG_FILE"
        
        python3 "${SCRIPTS_DIR}/create_variant_matrix_v2.py" \
            -s "$somatic_file" \
            -d "$donor_dir" \
            -o "${donor_name}_variant_sparse_matrix_withoutstrand.tsv.gz" \
            --workers "$WORKERS" 2>&1 | tee -a "$LOG_FILE"
        
        STEP4_EXIT=$?
    fi
    
    # 记录完成状态
    {
    echo -e "\n========================================"
    echo "完成处理: $donor_name"
    echo "完成时间: $(date)"
    echo "文件列表:"
    ls -lh *.tsv *.tsv.gz 2>/dev/null || echo "  无输出文件"
    echo "========================================"
    } | tee -a "$LOG_FILE"
    
done

# 最终汇总
{
echo -e "\n========================================"
echo "所有处理完成!"
echo "完成时间: $(date)"
echo "处理的Donor总数: ${#DONOR_DIRS[@]}"
echo "输出目录: $OUTPUT_DIR"
echo "日志文件: $LOG_FILE"
echo "各Donor输出:"
for donor_dir in "${DONOR_DIRS[@]}"; do
    donor_name=$(basename "$donor_dir")
    donor_output="${OUTPUT_DIR}/${donor_name}"
    if [ -d "$donor_output" ]; then
        file_count=$(find "$donor_output" -name "*.tsv" -o -name "*.tsv.gz" 2>/dev/null | wc -l)
        echo "  - $donor_name: $file_count 个文件"
    fi
done
echo "========================================"
} | tee -a "$LOG_FILE"

echo "处理完成! 详细日志见: $LOG_FILE"