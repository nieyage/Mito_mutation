# rawdata pwd: /data/R04/liyh526/project/00_PLOG_aging/01_data/00_raw_data/polg_27month/plog_mouse-BMMC_joint

# 1. QC for fastq

#!/bin/bash
#PBS -N fastqc_analysis
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=8G

cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger
ln -s /data/R04/liyh526/project/00_PLOG_aging/01_data/00_raw_data/polg_27month/plog_mouse-BMMC_joint/BMMC_27m_ATAC_H10 .

# 运行 FastQC 进行质量检测
mkdir -p /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/fastqc_results

cat BMMC_27m_ATAC_H10.input | while read id
do
    arr=($id)
    sample=${arr[0]}
    input1=${arr[1]}
    input2=${arr[2]}
    
    echo "Running FastQC for sample: $sample"
    fastqc ./BMMC_27m_ATAC_H10/*S[5-8]*R1*gz ./BMMC_27m_ATAC_H10/*S[5-8]*R2*gz ./BMMC_27m_ATAC_H10/*S[5-8]*R3*gz -o /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/fastqc_results -t 16
done


# only for S5-S8


# 检查 FastQC 结果，特别关注 100-150bp 区域的质量下降
#!/bin/bash
#PBS -N check_fastq_quality
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=8G

cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/

# 创建结果目录
mkdir -p quality_summary

echo "Analyzing quality drop in 100-150bp region..."
cat > check_quality.sh << 'EOF'
#!/bin/bash
sample=$1
r1_file="fastqc_results/${sample}_R1_fastqc/fastqc_data.txt"
r2_file="fastqc_results/${sample}_R2_fastqc/fastqc_data.txt"

check_quality_drop() {
    local file=$1
    local read_end=$2
    
    if [ ! -f "$file" ]; then
        echo "File not found: $file"
        return 1
    fi
    
    # 提取per base sequence quality部分
    awk '/>>Per base sequence quality/,/>>END/' "$file" | \
        grep -v "^>>" | grep -v "^#" | grep -v "^$" > /tmp/quality_data.txt
    
    # 计算前100bp的平均质量和100-150bp的平均质量
    pre_100_avg=$(awk '$1+0 <= 100 {sum+=$2; count++} END {print sum/count}' /tmp/quality_data.txt 2>/dev/null)
    post_100_avg=$(awk '$1+0 > 100 && $1+0 <= 150 {sum+=$2; count++} END {print sum/count}' /tmp/quality_data.txt 2>/dev/null)
    
    if [ -z "$pre_100_avg" ] || [ -z "$post_100_avg" ]; then
        echo "Cannot calculate quality for $read_end"
        return 1
    fi
    
    echo "  $read_end: 1-100bp avg quality: $pre_100_avg"
    echo "  $read_end: 101-150bp avg quality: $post_100_avg"
    
    # 判断是否有质量下降（下降超过2个质量分数）
    if (( $(echo "$pre_100_avg - $post_100_avg > 2" | bc -l) )); then
        echo "  WARNING: Quality drop detected in $read_end (difference: $(echo "$pre_100_avg - $post_100_avg" | bc))"
        return 0  # 返回0表示需要修剪
    else
        echo "  OK: No significant quality drop in $read_end"
        return 1  # 返回1表示不需要修剪
    fi
}

echo "=== Quality report for $sample ==="
needs_trim=0

if check_quality_drop "$r1_file" "R1"; then
    needs_trim=1
fi

if check_quality_drop "$r2_file" "R2"; then
    needs_trim=1
fi

if [ $needs_trim -eq 1 ]; then
    echo "RECOMMENDATION: Trim last 50bp for $sample"
else
    echo "RECOMMENDATION: No need to trim for $sample"
fi
echo ""
EOF

chmod +x check_quality.sh

# 3. 检查每个样本
echo "Generating quality report..."
cat rawdata/AR3_C5.input | while read id
do
    arr=($id)
    sample=${arr[0]}
    
    echo "Checking $sample..."
    ./check_quality.sh $sample >> quality_summary/quality_report.txt
done

# 4. 总结需要修剪的样本
echo "=== Summary ==="
echo "Samples needing trim (quality drop >2 in 100-150bp region):"
grep -B2 "RECOMMENDATION: Trim" quality_summary/quality_report.txt | grep "for " | awk '{print $NF}'

echo "Samples not needing trim:"
grep -B2 "RECOMMENDATION: No need" quality_summary/quality_report.txt | grep "for " | awk '{print $NF}'

echo "Detailed report saved to: quality_summary/quality_report.txt"





# 检查 multiqc 报告或手动查看 FastQC 结果
# 如果发现 100-150bp 质量下降，执行修剪步骤

# 判断标准：通常查看 per_base_sequence_quality 模块
# 如果 100-150bp 的质量分数明显低于前 100bp，建议修剪


# 2. trim the last 50 bp

#!/bin/bash
#PBS -N trim_fastq
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=10G

cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/

mkdir -p trim_data

cat BMMC_27m_ATAC_H10.input | while read id
do
    arr=($id)
    sample=${arr[0]}
    input1=${arr[1]}
    input2=${arr[2]}
    
    echo "Processing sample: $sample"
    
    # 使用 Trimmomatic 修剪最后50bp
    # CROP:100 保留前100bp，去掉后面的部分
    java -jar /md01/nieyg/ori/biosoft/package/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -phred33 \
        ./BMMC_27m_ATAC_H10/$input1 ./BMMC_27m_ATAC_H10/$input2 \
        ./trim_data/$input1 ./trim_data/$sample_unpaired_R1.fq.gz \
        ./trim_data/$input2 ./trim_data/$sample_unpaired_R3.fq.gz \
        CROP:100 LEADING:35 TRAILING:35 CROP:100 MINLEN:35 AVGQUAL:30 -threads 12
    
    echo "Finished trimming for $sample"
done

# check the fastqc of trimed data

cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/trim_data
fastqc *gz 
multiqc 




# 3. fastq_pair the R2(barcode) with R1,R3(reads)
#!/bin/bash
#PBS -N pair_fastq
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=8G

cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/trim_data

# 定义样本名称（根据您的实际情况修改）
sample_names=("S5" "S6" "S7" "S8")

for sample_name in "${sample_names[@]}"; do
    echo "Processing sample: $sample_name"
    
    # Decompress R1 and R2 FASTQ files
    gunzip -c "/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/BMMC_27m_ATAC_H10/BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq.gz" > "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq"
    gunzip -c "BMMC_27m_ATAC_H10_${sample_name}_L001_R1_001.fastq.gz" > "BMMC_27m_ATAC_H10_${sample_name}_L001_R1_001.fastq"
 
    # Use fastq_pair to pair the reads
    fastq_pair "BMMC_27m_ATAC_H10_${sample_name}_L001_R1_001.fastq" "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq"

    # Rename the paired file
    mv "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq.paired.fq" "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq"
    gzip "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq"

    # Remove temporary files
    rm "BMMC_27m_ATAC_H10_${sample_name}_L001_R2_001.fastq"
    rm "BMMC_27m_ATAC_H10_${sample_name}_L001_R1_001.fastq"
    
    echo "Finished pairing for $sample_name"
done


#!/bin/bash

# List of sample names to process
sample_names=("S5" "S6" "S7" "S8")

for sample_name in "${sample_names[@]}"; do

    # Decompress R1 and R2 FASTQ files
    gunzip -c "/data/R03/zhangwx/rawdata/mouse/Yage-mouse/AR3/AR3_C4/AR3_C4_${sample_name}_L001_R2_001.fastq.gz" > "AR3_C4_${sample_name}_L001_R2_001.fastq"
    gunzip -c "AR3_C4_${sample_name}_L001_R1_001.fastq.gz" > "AR3_C4_${sample_name}_L001_R1_001.fastq"

    # Use fastq_pair to pair the reads
    fastq_pair "AR3_C4_${sample_name}_L001_R1_001.fastq" "AR3_C4_${sample_name}_L001_R2_001.fastq"

    # Remove temporary files
    rm "AR3_C4_${sample_name}_L001_R2_001.fastq"
    rm "AR3_C4_${sample_name}_L001_R1_001.fastq"

    # Rename the paired file
    mv "AR3_C4_${sample_name}_L001_R2_001.fastq.paired.fq" "AR3_C4_${sample_name}_L001_R2_001.fastq"

    # Compress the paired file
    gzip "AR3_C4_${sample_name}_L001_R2_001.fastq"
done


# 2. Cellranger (STAR) by masked genome (trimed fastq)

# ref: /md01/nieyg/ref/hard-mask/mm10_hard_masked
# libraries.csv
fastqs,sample,library_type
/data/R04/liyh526/project/00_PLOG_aging/01_data/00_raw_data/polg_27month/plog_mouse-BMMC_joint/BMMC_27m_RNA_H3,BMMC_27m_RNA_H3,Gene Expression
/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/trim_data,BMMC_27m_ATAC_H10,Chromatin Accessibility


#PBS -N cellranger_masked
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger
cellranger-arc count --id=BMMC_27m_ATAC_H10 \
                       --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                       --libraries=/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/libraries.csv \
                       --localcores=24 \
                       --localmem=64
