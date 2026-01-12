# Step1: make the shifted ref 

sh make_shifted_ref.sh 

# ref:https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
bwa index mm10.chrM.fasta
bwa index mm10.chrM.shifted_by_8000_bases.fasta
samtools faidx Homo_sapiens_assembly38.chrM.fasta
samtools faidx Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=mm10.chrM.fasta O=mm10.chrM.dict
java -Xmx4g -jar /public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=mm10.chrM.shifted_by_8000_bases.fasta O=mm10.chrM.shifted_by_8000_bases.dict


# Step2: Get unmapped and chrM reads and Remapped with shifted chrM genome by bwa
input_bam=/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/BMMC_27m_ATAC/outs/atac_possorted_bam.bam
samtools view -b $input_bam chrM > BMMC_27m_chrM_mapped.bam
samtools view -b -f 4 $input_bam > BMMC_27m_unmapped.bam
samtools merge BMMC_27m_possorted_chrM_and_unmapped.bam BMMC_27m_unmapped.bam BMMC_27m_chrM_mapped.bam
samtools index BMMC_27m_possorted_chrM_and_unmapped.bam
rm BMMC_27m_unmapped.bam BMMC_27m_chrM_mapped.bam

# https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/mitochondria_references

unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.shifted_by_8000_bases.fasta
BMMC_27m_chrM_unmapped_bam=/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/BMMC_27m_possorted_chrM_and_unmapped.bam
# remapping by shifted genome 
samtools collate -Oun128 $BMMC_27m_chrM_unmapped_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $BMMC_27m_chrM_unmapped_bam|grep ^@RG) $shifted_chrM_ref - \
  | samtools sort -@4 -m4g -o BMMC_27m_chrM_unmapped_bwa_shifted.bam - 

samtools index BMMC_27m_chrM_unmapped_bwa_shifted.bam

#(sleep 3600 && ./remap_chrM_shift.sh) &

nohup sh -c './remap_chrM_shift.sh' > output.log 2>&1 &

nohup samtools sort -@ 12 -t CB BMMC_27m_chrM_unmapped_bwa_shifted.bam -o BMMC_27m_chrM_unmapped_bwa_shifted_sorted_CB.bam & 
nohup samtools sort -@ 36 -t CB BMMC_27m_possorted_chrM_and_unmapped.bam -o BMMC_27m_chrM_unmapped_bwa_unshifted_sorted_CB.bam & 

# R get high quality cells 

tmp<- readRDS("/data/R04/liyh526/project/00_PLOG_aging/01_data/02_processed_data/BMMC_PBMC34_clear_celltype.rds")
data<- tmp@meta.data
data<-data[which(data$orig.ident=="27M_BMMC"),]
anno_data<- data.frame(barcode=data$atac_barcode,celltype=data$cell_type)
write.csv(anno_data,"/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/barcode_celltype.csv",row.names=F)
anno_data<- data.frame(barcode=data$gex_barcode,celltype=data$cell_type)
write.csv(anno_data,"/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv",row.names=F)

cut -d, -f1 /md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv | tr -d '"' > barcodes.txt

# split bam by barcode 
nohup ./splitbam_unshifted.py -i BMMC_27m_chrM_unmapped_bwa_unshifted_sorted_CB.bam -b barcodes.txt > splitbam_unshifted_output.log 2>&1 &
nohup ./splitbam_shifted.py -i BMMC_27m_chrM_unmapped_bwa_shifted_sorted_CB.bam -b barcodes.txt > splitbam_shifted_output.log 2>&1 &

# for test


nohup sh ./SNV_Calling.sh > SNV_Calling_output3.log 2>&1 &
