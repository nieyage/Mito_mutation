# Step1: Get unmapped and chrM reads
# input_bam= /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
input_bam=/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam

samtools view -b $input_bam chrM > PBMC_chrM_mapped.bam
samtools view -b -f 4 $input_bam > PBMC_unmapped.bam
samtools merge PBMC_possorted_chrM_and_unmapped.bam PBMC_unmapped.bam PBMC_chrM_mapped.bam
samtools index PBMC_possorted_chrM_and_unmapped.bam
rm PBMC_unmapped.bam PBMC_chrM_mapped.bam

# Step2: Remapped with shifted/unshifted chrM genome by bwa
# ref:https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
bwa index Homo_sapiens_assembly38.chrM.fasta
bwa index Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
samtools faidx Homo_sapiens_assembly38.chrM.fasta
samtools faidx Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=Homo_sapiens_assembly38.chrM.fasta O=Homo_sapiens_assembly38.chrM.dict
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta O=Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict

# https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references
unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
PBMC_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_possorted_chrM_and_unmapped.bam

nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $unshifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_unshifted.bam - &

# 2. remapping by shifted genome 
nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_shifted.bam - &

samtools index PBMC_chrM_unmapped_bwa_unshifted.bam
samtools index PBMC_chrM_unmapped_bwa_shifted.bam

# Step3: filter bam file 

samtools view -@ 12 -h -b -q 30 -f 2 -F 3852 PBMC_chrM_unmapped_bwa_unshifted.bam  -o PBMC_chrM_unmapped_bwa_unshifted_filtered1.bam
samtools index PBMC_chrM_unmapped_bwa_unshifted_filtered1.bam
python  para_filter_step2_V5.py \
    --bam PBMC_chrM_unmapped_bwa_unshifted_filtered1.bam \
    --out  PBMC_chrM_unmapped_bwa_unshifted_filtered2.bam \
    --region chrM:1-16569 \
    --threads 6 \
    --max_edits 3 \
    --min_reads_per_barcode 1000

samtools view -@ 12 -h -b -q 30 -f 2 -F 3852 PBMC_chrM_unmapped_bwa_shifted.bam  -o PBMC_chrM_unmapped_bwa_shifted_filtered1.bam
samtools index PBMC_chrM_unmapped_bwa_shifted_filtered1.bam

python  para_filter_step2_V5.py \
    --bam PBMC_chrM_unmapped_bwa_shifted_filtered1.bam \
    --out  PBMC_chrM_unmapped_bwa_shifted_filtered2.bam \
    --region chrM:1-16569 \
    --threads 6 \
    --max_edits 3 \
    --min_reads_per_barcode 1000

# Step4: change the position of shifted bam

#!/usr/bin/env python3
import pysam
import sys

def convert_bam_coordinates(shift_bam, unshift_bam, shift=8000, mt_length=16569, mt_chrom="chrM"):
    """
    将shifted BAM文件转换回原始坐标
    """
    # 打开输入BAM文件
    with pysam.AlignmentFile(shift_bam, "rb") as infile:
        # 创建输出BAM文件
        with pysam.AlignmentFile(unshift_bam, "wb", header=infile.header) as outfile:
            for read in infile:
                # 只处理线粒体染色体的read
                if read.reference_name == mt_chrom and not read.is_unmapped:
                    # 转换坐标
                    old_pos = read.reference_start + 1  # 转为1-based
                    if old_pos <= (mt_length - shift):
                        new_pos = old_pos + shift
                    else:
                        new_pos = old_pos - (mt_length - shift)
                    
                    # 设置新的位置（转回0-based）
                    read.reference_start = new_pos - 1
                    
                    # 如果有mate信息也需要转换
                    if read.is_paired and read.next_reference_name == mt_chrom and not read.mate_is_unmapped:
                        old_mpos = read.next_reference_start + 1
                        if old_mpos <= (mt_length - shift):
                            new_mpos = old_mpos + shift
                        else:
                            new_mpos = old_mpos - (mt_length - shift)
                        read.next_reference_start = new_mpos - 1
                
                outfile.write(read)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_bam.py shift_bam.bam unshift_bam.bam")
        sys.exit(1)
    
    shift_bam = sys.argv[1]
    unshift_bam = sys.argv[2]
    convert_bam_coordinates(shift_bam, unshift_bam)


python convert_bam.py PBMC_chrM_unmapped_bwa_shifted_filtered2.bam PBMC_chrM_unmapped_bwa_shifted_filtered2_convert_pos.bam
samtools sort -@ 12 PBMC_chrM_unmapped_bwa_shifted_filtered2_convert_pos.bam -o PBMC_chrM_unmapped_bwa_shifted_filtered2_convert_pos_sorted.bam
samtools index -@ 12 PBMC_chrM_unmapped_bwa_shifted_filtered2_convert_pos_sorted.bam

samtools sort -@ 12 PBMC_chrM_unmapped_bwa_unshifted_filtered2.bam -o PBMC_chrM_unmapped_bwa_unshifted_filtered2_sorted.bam
samtools index -@ 12 PBMC_chrM_unmapped_bwa_unshifted_filtered2_sorted.bam


# Step5: Merge the 2 bam file(shifted bam for Dloop and unshifted bam for other region)
shifted_bam_sorted=PBMC_chrM_unmapped_bwa_shifted_filtered2_convert_pos_sorted.bam
unshifted_bam_sorted=PBMC_chrM_unmapped_bwa_unshifted_filtered2_sorted.bam
# 从shifted_bam_sorted提取 1-575bp 和 16025-16569bp
samtools view -b $shifted_bam_sorted chrM:1-575 chrM:16025-16569 > part1.bam

# 去除起始位置在576-16024之间的reads
samtools view -h part1.bam | awk '
    /^@/ {print; next}  # 保留头信息
    $4 < 576 || $4 >= 16025 {print}  # 只输出起始位置<576或>=16025的reads
' | samtools view -b -o filtered_part1.bam

# samtools view filtered_part1.bam | awk '{print $4}' | sort -n | head -5
# samtools view filtered_part1.bam | awk '{print $4}' | sort -n | tail -5

# 从unshifted_bam_sorted提取 576-16024bp
samtools view -@ 12 -b $unshifted_bam_sorted chrM:576-16024 > part2.bam
samtools view -@ 12 -h part2.bam | awk '
    /^@/ {print; next}  # 保留头信息
    $4 >= 576 {print}   # 只输出起始位置 >= 576 的reads
' | samtools view -@ 12 -b -o part2_filtered.bam

# samtools view part2_filtered.bam | awk '{print $4}' | sort -n | head -5
# samtools view part2_filtered.bam | awk '{print $4}' | sort -n | tail -5

samtools merge -@ 12 PBMC_chrM_unmapped_bwa_merged.bam filtered_part1.bam part2_filtered.bam
samtools sort -@ 12 PBMC_chrM_unmapped_bwa_merged.bam -o PBMC_chrM_unmapped_bwa_merged_sorted.bam
samtools index -@ 12 PBMC_chrM_unmapped_bwa_merged_sorted.bam

rm part1.bam part2.bam filtered_part1.bam part2_filtered.bam PBMC_chrM_unmapped_bwa_merged.bam
rm -rf filter_bam_temp 
rm *filtered1.bam* 
rm *filtered2.bam* 

# Step6: count the matrix (h5)
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/02_getchrM_unmapped_masked/PBMC_chrM_unmapped_bwa_merged_sorted.bam .
python sumstatsBP_parallel_hdf5_workers.py PBMC_chrM_unmapped_bwa_merged_sorted.bam  output.h5 16569 20 30 True 8 temp 1>log 2>err


library(systemfonts)


