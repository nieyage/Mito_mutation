# Step1: Get unmapped and chrM reads
# input_bam= /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
input_bam=/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam
samtools view -b $input_bam chrM > PBMC_chrM_mapped.bam
samtools view -b -f 4 $input_bam > PBMC_unmapped.bam
samtools merge PBMC_possorted_chrM_and_unmapped.bam PBMC_unmapped.bam PBMC_chrM_mapped.bam
samtools index PBMC_possorted_chrM_and_unmapped.bam
rm PBMC_unmapped.bam 

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
PBMC_chrM_unmapped_bam=PBMC_possorted_chrM_and_unmapped.bam
nohup samtools collate -Oun128 $PBMC_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $PBMC_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o PBMC_chrM_unmapped_bwa_shifted.bam - &

# filtered unmapped reads 
samtools index PBMC_chrM_unmapped_bwa_shifted.bam
samtools view -@ 12 -b PBMC_chrM_unmapped_bwa_shifted.bam chrM > PBMC_chrM_unmapped_bwa_shifted2.bam

# Step3: filter bam 
shift_bam=/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/PBMC_chrM_unmapped_bwa_shifted2.bam
unshift_bam=/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/PBMC_chrM_mapped.bam

# unshift_bam
samtools view -@ 12 -h -b -q 30 -f 2 -F 3852 $unshift_bam  -o PBMC_chrM_mapped_filtered1.bam
samtools index PBMC_chrM_mapped_filtered1.bam
python  para_filter_step2_V5.py \
    --bam PBMC_chrM_mapped_filtered1.bam \
    --out  PBMC_chrM_mapped_filtered2.bam \
    --region chrM:1-16569 \
    --threads 6 \
    --max_edits 3 \
    --min_reads_per_barcode 1000

# shift_bam
samtools view -@ 12 -h -b -q 30 -f 2 -F 3852 $shift_bam  -o PBMC_chrM_shift_filtered1.bam
samtools index PBMC_chrM_shift_filtered1.bam
python  para_filter_step2_V5.py \
    --bam PBMC_chrM_shift_filtered1.bam \
    --out  PBMC_chrM_shift_filtered2.bam \
    --region chrM:1-16569 \
    --threads 6 \
    --max_edits 3 \
    --min_reads_per_barcode 1000
# only keep the DLoop region chrM:7000-10000 
samtools index PBMC_chrM_shift_filtered2.bam 
samtools view -h -b PBMC_chrM_shift_filtered2.bam chrM:7000-10000 > PBMC_chrM_Dloop_filtered2.bam

# Step4: sort by CB 
all_bam=/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/PBMC_chrM_mapped_filtered2.bam
Dloop_bam=/md01/nieyg/project/mito_mutation/01_pipeline/08_v4/PBMC_chrM_Dloop_filtered2.bam

samtools sort -@ 36 -t CB $all_bam -o PBMC_chrM_sorted_CB.bam
samtools sort -@ 36 -t CB $Dloop_bam -o PBMC_chrM_Dloop_sorted_CB.bam

# input a barcode,celltype annotation file 
cut -d, -f1 /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv | tr -d '"' > barcodes.txt

# Step5: split bam by barcode 
nohup python ./splitbam_unshifted.py -i PBMC_chrM_sorted_CB.bam -b barcodes.txt > splitbam_unshifted_output.log 2>&1 &
nohup python ./splitbam_shifted.py -i PBMC_chrM_Dloop_sorted_CB.bam -b barcodes.txt > splitbam_shifted_output.log 2>&1 &

# Step6: mutation calling for each barcode;
nohup sh ./SNV_calling_v1.sh > SNV_Calling_outputv1.log 2>&1 &



nohup sh ./SNV_calling_v3.sh > SNV_Calling_output.log 2>&1 &



nohup sh ./SNV_calling_v4.sh > SNV_Calling_output_v4.log 2>&1 &
