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

# 检查 multiqc 报告或手动查看 FastQC 结果
# 如果发现 100-150bp 质量下降，执行修剪步骤
# 判断标准：通常查看 per_base_sequence_quality 模块, 如果 100-150bp 的质量分数明显低于前 100bp，建议修剪

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

# 定义样本名称
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


# 4. Cellranger (STAR) by masked genome (trimed fastq)

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
