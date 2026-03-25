# Step1: Get unmapped and chrM reads
# input_bam= /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
input_bam= /md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/atac_possorted_bam.bam

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

# 1. remapping by unshifted genome 

unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
PBMC_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_possorted_chrM_and_unmapped.bam

nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $unshifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_unshifted.bam - &

# 2. remapping by shifted genome 
nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_shifted.bam - &

samtools index PBMC_chrM_unmapped_bwa_unshifted.bam
samtools index PBMC_chrM_unmapped_bwa_shifted.bam

# Step3: Realign by GATK

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R $unshifted_chrM_ref  \
-T RealignerTargetCreator  -nt 32 \
-I PBMC_chrM_unmapped_bwa_unshifted.bam \
-o PBMC_chrM_unmapped_bwa_unshifted.target.intervals

nohup java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R $unshifted_chrM_ref  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I PBMC_chrM_unmapped_bwa_unshifted.bam \
-targetIntervals PBMC_chrM_unmapped_bwa_unshifted.target.intervals \
-o PBMC_chrM_unmapped_bwa_unshifted.realign.bam &

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R $shifted_chrM_ref  \
-T RealignerTargetCreator  -nt 32 \
-I PBMC_chrM_unmapped_bwa_shifted.bam \
-o PBMC_chrM_unmapped_bwa_shifted.target.intervals

nohup java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R $shifted_chrM_ref  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I PBMC_chrM_unmapped_bwa_shifted.bam \
-targetIntervals PBMC_chrM_unmapped_bwa_shifted.target.intervals \
-o PBMC_chrM_unmapped_bwa_shifted.realign.bam &

nohup samtools sort -t CR PBMC_chrM_unmapped_bwa_unshifted.realign.bam -o PBMC_chrM_unmapped_bwa_unshifted_realign_sorted.bam &
nohup samtools sort -t CR PBMC_chrM_unmapped_bwa_shifted.realign.bam -o PBMC_chrM_unmapped_bwa_shifted_realign_sorted.bam &


cd ../03_split_bam
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_unshifted_realign_sorted.bam .
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_shifted_realign_sorted.bam .
python split_bam_shift_unshift.py





