# get
input_fasta=/md01/nieyg/ref/hard-mask/genome_modify_hg38/modify_genome.fa
# 直接从文件中提取 chrM（支持压缩文件）
seqkit grep -p "chrM" $input_fasta > chrM.fa
# 提取环状区域 16519-150，并保存

# dloop_flank_16569_150bp.fa

ACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACAT
CTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGAC
ATCACGATG

#dloop_flank_1_150bp.fa

GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATC

# merge
>dloop_flank
ACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACAT
CTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGAC
ATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATT
TGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCA
CCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATC


cat /md01/nieyg/ref/hard-mask/genome_modify_hg38/modify_genome.fa dloop_flank.fa > add_flank.fasta



input_gtf=/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf
input_fasta=/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_input/add_flank.fasta


# 2) create cellranger genome 
# hg38_hard_masked_addflank.config
organism: "human" # same as before
genome: ["hg38_hard_masked_addflank"] # same as before
input_fasta: ["/md01/nieyg/ref/hard-mask/dnbc4tools_input/hg38_input/add_flank.fasta"] # updated from bedtools maskfasta with masked regions and add flank
input_gtf: ["/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"] # same as before
#non_nuclear_contigs: ["chrM"] # same as before
#input_motifs: "/md01/nieyg/ref/hard-mask/genome_modify/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt" # same as before

nohup cellranger-arc mkref --config=hg38_hard_masked_addflank.config &
