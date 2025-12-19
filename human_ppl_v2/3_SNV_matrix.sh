# Step1: count the matrix (h5)
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_chrM_unmapped_bwa_merged_sorted.bam .
ln -s /md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped_masked/PBMC_chrM_unmapped_bwa_merged_sorted.bam.bai .

python sumstatsBP_parallel_hdf5_workers.py PBMC_chrM_unmapped_bwa_merged_sorted.bam  output.h5 16569 20 30 True 8 temp 1>log 2>err


unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta

python /md01/jinxu/Project/mgatk-speedup/22_VarCall/2_variantCall.py  output.h5 PBMC 16569 1 /md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta

/md01/jinxu/Project/mgatk-speedup/

