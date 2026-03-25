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
python SNV_summary_v2.py -i /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_allcell -o ./snv.tsv -v 


# 2) get germline or somatic snv

python3 filter_mutations.py -h

python3 filter_mutations.py \
  -i snv.tsv.gz \
  -o ./ \
  --germline-output germline_withoutstrand.tsv \
  --somatic-output somatic_withoutstrand.tsv \
  --strand-min 0 \
  --strand-max 1 \
  --min-depth 10 \
  --alt-ratio-threshold 0.90 \
  --no-compress


# 3. variant_sparse_matrix.tsv

python3 create_variant_matrix_v2.py -s somatic_withoutstrand.tsv \
-d /md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_allcell \
-o variant_sparse_matrix_withoutstrand.tsv.gz --workers 10


# For mix sample 
# split by donor 


python organize_files_parallel.py \
  -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib1_output/SNV_Calling_allcell \
  -a /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation_new/lib1_annotation.csv \
  -o /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib1_output/output \
  -w 16 

python organize_files_parallel.py \
  -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib3_output/SNV_Calling_allcell \
  -a /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation_new/lib3_annotation.csv \
  -o /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib3_output/output \
  -w 36  

python organize_files_parallel.py \
  -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/SNV_Calling_allcell \
  -a /data/R01/renwx5/MT_mution/mutation_outs/sample_annotation_new/lib5_annotation.csv \
  -o /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/output \
  -w 16  # 使用16个工作线程

python organize_files_parallel.py \
  -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/SNV_Calling_allcell \
  -a /md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv \
  -o /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/output \
  -w 36  # 使用16个工作线程

bash /md01/nieyg/pipeline/mito_mutation/01_human/process_donors.sh \
  --base-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib1_output/output \
  --output-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib1_output/donor_results \
  --threads 32 \
  --workers 20

bash /md01/nieyg/pipeline/mito_mutation/01_human/process_donors.sh \
  --base-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib3_output/output \
  --output-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib3_output/donor_results \
  --threads 32 \
  --workers 20


bash /md01/nieyg/pipeline/mito_mutation/01_human/process_donors.sh \
  --base-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/output \
  --output-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results \
  --threads 32 \
  --workers 20

bash /md01/nieyg/pipeline/mito_mutation/01_human/process_donors.sh \
  --base-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/output \
  --output-dir /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/donor_results \
  --threads 32 \
  --workers 20


