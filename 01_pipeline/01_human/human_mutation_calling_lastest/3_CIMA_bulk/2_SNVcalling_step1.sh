# female_57 

annotation_file=/data/R01/renwx5/MT_mution/CIMA/annotation/CNR1281138_annotation.csv
all_bam=/data/R01/renwx5/MT_mution/CIMA/align/new_ref/CNR1281138/outs/alignment.fragments.sorted.tagged.bam

# Step1: Get unmapped and chrM reads
input_bam=$all_bam
samtools view -b $input_bam chrM > chrM_mapped.bam
samtools view -b $input_bam dloop_flank > dloop_flank_mapped.bam
samtools merge chrM_and_dloop_flank.bam dloop_flank_mapped.bam chrM_mapped.bam
samtools index chrM_and_dloop_flank.bam

cut -d',' -f1 $annotation_file > barcode.txt

subset_bam_tool="/md01/nieyg/software/subset-bam_linux"
/md01/nieyg/software/subset-bam_linux --bam chrM_and_dloop_flank.bam \
        --cell-barcodes ./barcode.txt \
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

samtools index chrM_unshifted.bam
samtools index Dloop_shifted.bam

# Step3:  split bam by celltype 
all_bam=/md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/chrM_unshifted.bam
Dloop_bam=/md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/Dloop_shifted.bam

mkdir celltype_barcode && cd celltype_barcode
awk -F',' 'NR>1 {print $1 >> $3".barcodes.txt"}' $annotation_file

mkdir unshifted_bam
mkdir Dloop_bam

for bc in *.barcodes.txt; do
    celltype=${bc%.barcodes.txt}
    /md01/nieyg/software/subset-bam_linux --bam $all_bam --cell-barcodes "$bc" --cores 16 --out-bam "./unshifted_bam/${celltype}.bam"
done


for bc in *.barcodes.txt; do
    celltype=${bc%.barcodes.txt}
    /md01/nieyg/software/subset-bam_linux --bam $Dloop_bam --cell-barcodes "$bc" --cores 16 --out-bam "./Dloop_bam/${celltype}.bam"
done


# Step4: mutation calling for each barcode;
# SNV_calling.sh中参数需要修改
nohup sh ./SNV_calling_bulk.sh > SNV_Calling_output_allcell.log 2>&1 &







