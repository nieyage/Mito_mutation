# input a barcode with one bead
awk -F',' '$8==1 && $9 !~ /;/ {print $1}' /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/outs/singlecell.csv > barcode.txt

# Step1: Get unmapped and chrM reads
input_bam=/md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/outs/alignment.fragments.sorted.tagged.bam
samtools view -b $input_bam chrM > chrM_mapped.bam
samtools view -b $input_bam dloop_flank > dloop_flank_mapped.bam

# samtools view -b -f 4 $input_bam > PBMC_unmapped.bam
samtools merge chrM_and_dloop_flank.bam dloop_flank_mapped.bam chrM_mapped.bam
samtools index chrM_and_dloop_flank.bam

subset_bam_tool="/md01/nieyg/software/subset-bam_linux"
/md01/nieyg/software/subset-bam_linux --bam chrM_and_dloop_flank.bam \
        --cell-barcodes /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/outs/barcode.txt \
        --cores 16 \
        --out-bam chrM_and_dloop_flank_filteredBC.bam
    


# rm chrM_mapped.bam dloop_flank_mapped.bam 

# Step2: Remapped with unshifted/shifted chrM genome by bwa
chrM_flank_bam=chrM_and_dloop_flank_filteredBC.bam

#  remapping by unshifted genome 
unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
nohup samtools collate -Oun128 $chrM_flank_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $chrM_flank_bam|grep ^@RG) $unshifted_chrM_ref - \
  | samtools sort -@4 -m4g -o chrM_flank_bwa_unshifted.bam - &
# filtered unmapped reads 
samtools index chrM_flank_bwa_unshifted.bam
samtools view -@ 12 -b chrM_flank_bwa_unshifted.bam chrM > chrM_unshifted.bam

#  remapping by shifted genome 
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
nohup samtools collate -Oun128 $chrM_flank_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $chrM_flank_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o chrM_flank_bwa_shifted.bam - &
# filtered unmapped reads 
samtools index chrM_flank_bwa_shifted.bam
samtools view -@ 12 -b chrM_flank_bwa_shifted.bam chrM > chrM_flank_bwa_shifted2.bam
# only keep the DLoop region chrM:7000-10000 
samtools index chrM_flank_bwa_shifted.bam
samtools view -h -b chrM_flank_bwa_shifted.bam chrM:7000-10000 > Dloop_shifted.bam

# Step3:  sort by CB 
all_bam=chrM_unshifted.bam
Dloop_bam=Dloop_shifted.bam
samtools sort -@ 36 -t CB $all_bam -o chrM_unshifted_sorted_CB.bam
samtools sort -@ 36 -t CB $Dloop_bam -o Dloop_shifted_sorted_CB.bam


# Step5: split bam by barcode 
python /md01/nieyg/pipeline/mito_mutation/01_human/split_bam.py -i chrM_unshifted_sorted_CB.bam -b /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/outs/barcode.txt -o ./unshifted_bam/
python /md01/nieyg/pipeline/mito_mutation/01_human/split_bam.py -i Dloop_shifted_sorted_CB.bam -b /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/outs/barcode.txt -o ./Dloop_bam/

# Step6: mutation calling for each barcode;
# SNV_calling.sh中参数需要修改
nohup sh ./SNV_calling.sh > SNV_Calling_output_allcell.log 2>&1 &



