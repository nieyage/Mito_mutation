# Step1: Get unmapped and chrM reads
# input_bam= /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
input_bam=/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam

samtools view -b $input_bam chrM > PBMC_chrM_mapped.bam
samtools view -b -f 4 $input_bam > PBMC_unmapped.bam
samtools merge PBMC_possorted_chrM_and_unmapped.bam PBMC_unmapped.bam PBMC_chrM_mapped.bam
samtools index PBMC_possorted_chrM_and_unmapped.bam
rm PBMC_unmapped.bam PBMC_chrM_mapped.bam


# Step2: Remapped with shifted chrM genome by bwa
# ref:https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
bwa index Homo_sapiens_assembly38.chrM.fasta
bwa index Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
samtools faidx Homo_sapiens_assembly38.chrM.fasta
samtools faidx Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta

# https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references

#  remapping by shifted genome 

unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
PBMC_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_possorted_chrM_and_unmapped.bam

nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_shifted.bam - &

samtools index PBMC_chrM_unmapped_bwa_shifted.bam

cd /md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3 
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_possorted_chrM_and_unmapped.bam .
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_chrM_unmapped_bwa_shifted.bam .
# Step3: sort by CB 
nohup samtools sort -@ 36 -t CB PBMC_chrM_unmapped_bwa_shifted.bam -o PBMC_chrM_unmapped_bwa_shifted_sorted_CB.bam & 
nohup samtools sort -@ 36 -t CB PBMC_possorted_chrM_and_unmapped.bam -o PBMC_chrM_unmapped_bwa_unshifted_sorted_CB.bam & 

cut -d, -f1 /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv | tr -d '"' > barcodes.txt

# Step4: split bam by barcode 
nohup ./splitbam_unshifted.py -i PBMC_chrM_unmapped_bwa_unshifted_sorted_CB.bam -b barcodes.txt > splitbam_unshifted_output.log 2>&1 &
nohup ./splitbam_shifted.py -i PBMC_chrM_unmapped_bwa_shifted_sorted_CB.bam -b barcodes.txt > splitbam_shifted_output.log 2>&1 &

# Step5: mutation calling for each barcode;

nohup sh ./SNV_Calling.sh > SNV_Calling_output3.log 2>&1 &



