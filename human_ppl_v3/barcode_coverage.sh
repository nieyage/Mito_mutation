#!/bin/bash

# 设置路径
input_dir="/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell"
output_file="barcode_coverage.tsv.gz"

echo "开始从counts文件中提取覆盖度信息..."
echo "使用awk批量处理..."
echo -e "Barcode\tPosition\tCount" | gzip > "${output_file}"

for counts_file in "${input_dir}"/*.counts; do
    if [[ ! -f "${counts_file}" ]]; then
        continue
    fi
    
    # 提取barcode名称（去掉.counts后缀）
    barcode=$(basename "${counts_file}" .counts)
    
    echo "处理: ${barcode}"
    
    # 提取数据并添加到压缩文件中
    awk -v barcode="${barcode}" '
        BEGIN {OFS="\t"}
        /^chrM/ || $1 ~ /^[0-9XY]+$/ {
            # 第2列是position，第4列是Count
            print barcode, $2, $4
        }
    ' "${counts_file}" | gzip >> "${output_file}"
done

echo "提取完成！输出文件: ${output_file}"

exit 0

