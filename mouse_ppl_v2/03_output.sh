# 1. barcode_coverage.tsv.gz

# Barcode Position        Count
# CATAAGCTCAAGCCTG-1      1       7
# TATTTGGAGGCTACTG-1      1       9
# CGTGACATCCTTAGGG-1      1       7
python extract_coverage.py -h

# python ../extract_coverage.py -i /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_test_v5 -o test_barcode_coverage.tsv.gz -p 16

python extract_coverage.py -i /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_allcell -o barcode_coverage.tsv.gz -p 16 -b 1000



# 2. snv.somatic.tsv
# 1) summary all snv

python SNV_summary_v2.py --help

# 基本用法
python SNV_summary_v2.py -i /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_test_v5 -o ./snv.tsv -v 


# 2) get germline or somatic snv

python3 filter_mutations.py -h

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

python3 create_variant_matrix_v2.py -s somatic.tsv \
-d /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_test_v5 \
-o variant_sparse_matrix.tsv.gz --workers 10


