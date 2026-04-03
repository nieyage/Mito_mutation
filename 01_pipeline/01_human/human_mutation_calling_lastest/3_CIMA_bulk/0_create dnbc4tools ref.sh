# Doc:https://mgi-tech-bioinformatics.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/README.html

#### Download using `wget`
wget -O dnbc4tools-3.0.tar.gz "ftp://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.0.tar.gz"

tar -xzvf dnbc4tools-3.0.tar.gz

# Test basic functionality
./dnbc4tools --help
./dnbc4tools --version

# Test specific modules
./dnbc4tools rna --help
./dnbc4tools atac --help
./dnbc4tools vdj --help

dnbc4tools=/md01/nieyg/software/dnbc4tools3.0/dnbc4tools
input_gtf=/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf
input_fasta=/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_input/add_flank.fasta

$dnbc4tools tools mkgtf --ingtf $input_gtf --output human_genes.filter.gtf --type gene_type
filtered_input_gtf=/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_input/human_genes.filter.gtf

$dnbc4tools atac mkref \
  --fasta $input_fasta \
  --ingtf $filtered_input_gtf \
  --species hg38_masked_ATAC_addflank

# /md01/jinxu/Project/CNP6227
cd /md01/nieyg/project/mito_mutation/05_publish_dnbc4tools
dnbc4tools=/md01/nieyg/software/dnbc4tools3.0/dnbc4tools

nohup $dnbc4tools atac run \
  --name CNR1280203 \
  --fastqs /md01/jinxu/Project/CNP6227/CNR1280203/ \
  --genomeDir /md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_masked_ATAC_addflank \
  --threads 32 --need_bam & 

# 27G 样本 32 线程 19:02 开始 23:42:24 结束
# Analysis Finished Elapsed Time: 4:31:04



