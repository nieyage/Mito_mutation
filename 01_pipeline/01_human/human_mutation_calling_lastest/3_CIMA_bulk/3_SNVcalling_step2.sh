# 1. merge the output (snv and count)

python /data/R02/nieyg/pipeline/mito_mutation/01_human/counts2h5ad_v2.py \
-i /md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/SNV_calling/ \
-o ../All_cell_combined_counts_v2.h5ad

python /data/R02/nieyg/pipeline/mito_mutation/01_human/merge_snvs_v2.py \
-i /md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/SNV_calling/ \
-o ../All_cells_merged.snv.tsv


# 2. barcode_coverage.tsv.gz

python /data/R02/nieyg/pipeline/mito_mutation/01_human/extract_coverage.py \
-i /md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/celltype_barcode/SNV_calling/ \
-o barcode_coverage.tsv.gz -p 7 -b 1000

# 3. filter mutation

python /md01/nieyg/project/mito_mutation/06_bulk_CIMA/01_test/filter_mito.py -i All_cells_merged.snv.tsv -o celltype_filtered_results.csv

# 4. variant_sparse_matrix.tsv

python3 /data/R02/nieyg/pipeline/mito_mutation/01_human/extract_variant_matrix_from_h5ad.py \
  -s celltype_filtered_results.csv \
  -a All_cell_combined_counts_v2.h5ad \
  -o variant_sparse_matrix.tsv.gz


