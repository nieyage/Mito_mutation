# 1. merge the output (snv and count)

python counts2h5ad_v2.py -i /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/mutation/SNVcalling_allcell/ \
-o ../All_cell_combined_counts_v2.h5ad

python3 /data/R02/nieyg/pipeline/mito_mutation/01_human/merge_snvs_v2.py -i /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/mutation/SNVcalling_allcell/ -o ../All_cells_merged.snv.tsv


# 2. barcode_coverage.tsv.gz

python /data/R02/nieyg/pipeline/mito_mutation/01_human/extract_coverage.py \
-i /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/mutation/SNVcalling_allcell/ \
 -o barcode_coverage.tsv.gz -p 16 -b 1000


# 3. SNV summary

python3 /data/R02/nieyg/pipeline/mito_mutation/01_human/mt_mutation_stats.py \
  -s ./All_cells_merged.snv.tsv \
  -a ./All_cell_combined_counts_v2.h5ad \
  -o final_mutation_report.tsv

# python /data/R02/nieyg/pipeline/mito_mutation/01_human/SNV_summary_v2.py -i /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools/CNR1280832/mutation/SNVcalling_allcell/ -o ./snv.tsv -v 


# 4. get germline or somatic snv


python3 /data/R02/nieyg/pipeline/mito_mutation/01_human/filter_mt_mutations.py \
  -i final_mutation_report.tsv \
  --min-conf-cells 1 \
  --strand-min 0.3 \
  --strand-max 0.7 \
  --germline-vaf 0.9 \
  --somatic-vaf-max 0.9 \
  --output-prefix Sample_test

scp -r NCM/ nieyg@202.116.105.23:/public/home/nieyg/done_project


# python3 filter_mutations.py \
#   -i snv.tsv.gz \
#   -o ./ \
#   --germline-output germline_withoutstrand.tsv \
#   --somatic-output somatic_withoutstrand.tsv \
#   --strand-min 0 \
#   --strand-max 1 \
#   --min-depth 10 \
#   --alt-ratio-threshold 0.90 \
#   --no-compress

# 5. variant_sparse_matrix.tsv

python3 /data/R02/nieyg/pipeline/mito_mutation/01_human/extract_variant_matrix_from_h5ad.py \
  -s Sample_test_somatic.tsv \
  -a All_cell_combined_counts_v2.h5ad \
  -o variant_sparse_matrix.tsv.gz

# python3 create_variant_matrix_v2.py -s somatic_withoutstrand.tsv \
# -d /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_allcell \
# -o variant_sparse_matrix_withoutstrand.tsv.gz --workers 10

