# 1. barcode_coverage.tsv.gz

# Barcode Position        Count
# CATAAGCTCAAGCCTG-1      1       7
# TATTTGGAGGCTACTG-1      1       9
# CGTGACATCCTTAGGG-1      1       7
python extract_coverage.py -h

python extract_coverage.py -i /path/to/counts -o output.tsv.gz -p 16 -b 1000



# 2. snv.somatic.tsv

python3 snv_summary.py

python3 filter_mutations.py -h

python3 filter_mutations.py \
  -i snv_filtered.tsv.gz \
  -o ./ \
  --germline-output germline_filtered.tsv \
  --somatic-output somatic_filtered.tsv \
  --strand-min 0.3 \
  --strand-max 0.7 \
  --min-depth 10 \
  --alt-ratio-threshold 0.90 \
  --no-compress

python3 filter_mutations.py \
  -i snv.tsv.gz \
  -o ./ \
  --germline-output germline.tsv \
  --somatic-output somatic.tsv \
  --strand-min 0.3 \
  --strand-max 0.7 \
  --min-depth 10 \
  --alt-ratio-threshold 0.90 \
  --no-compress


# 3. variant_sparse_matrix.tsv
# python3 create_variant_matrix.py -s somatic_filtered.tsv \
# -d /md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell/ \
# -o variant_sparse_matrix_filtered.tsv.gz
# 
# python3 create_variant_matrix.py -s somatic.tsv \
# -d /md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell/ \
# -o variant_sparse_matrix.tsv.gz

python3 create_variant_matrix_v2.py -s somatic_filtered.tsv \
-d /md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell/ \
-o variant_sparse_matrix_filtered_v2.tsv.gz --workers 10

python3 create_variant_matrix_v2.py -s somatic.tsv \
-d /md01/nieyg/project/mito_mutation/01_pipeline/08_v4/masked_SNVcalling_percell_allcell/ \
-o variant_sparse_matrix_v2.tsv.gz --workers 10

# 4. plot 

Rscript mito_analysis.R \
  --coverage barcode_coverage.tsv.gz \
  --variants variant_sparse_matrix_filtered_v2.tsv.gz \
  --somatic somatic_filtered.tsv \
  --celltype /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv \
  --bam /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam \
  --barcodes /md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/barcodes.txt \
  --reference /md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta \
  --output-dir ./filtered_v2/ \
  --prefix filtered_v2

Rscript mito_analysis.R --help

Rscript mito_analysis.R \
  --coverage barcode_coverage.tsv.gz \
  --variants variant_sparse_matrix_v2.tsv.gz \
  --somatic somatic_filtered.tsv \
  --celltype /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv \
  --bam /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam \
  --barcodes /md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/barcodes.txt \
  --reference /md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta \
  --output-dir ./withoutfiltered_v2/ \
  --prefix withoutfiltered_v2



