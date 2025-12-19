# 1. create masked genome 

# blacklist from: https://github.com/caleblareau/mitoblacklist/blob/master/combinedBlacklist/hg38.full.blacklist.bed

# follow https://github.com/caleblareau/mgatk/wiki/Increasing-coverage-from-10x-processing

# pwd: /md01/nieyg/ref/mitoblacklist/mitoblacklist-master/combinedBlacklist/hg38.full.blacklist.bed

# 1) modify genome
mv genome.fa old_genome.fa
# mv genome.fa.flat old_genome.fa.flat
# mv genome.fa.gdx old_genome.fa.gdx
bedtools maskfasta -fi old_genome.fa -bed /md01/nieyg/ref/mitoblacklist/mitoblacklist-master/combinedBlacklist/hg38.full.blacklist.bed  -fo modify_genome.fa


# 2) create cellranger genome 
# hg38_hard_masked.config
organism: "human" # same as before
genome: ["hg38_hard_masked"] # same as before
input_fasta: ["/md01/nieyg/ref/hard-mask/genome_modify_hg38/modify_genome.fa"] # updated from bedtools maskfasta with masked regions
input_gtf: ["/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz"] # same as before
non_nuclear_contigs: ["chrM"] # same as before
#input_motifs: "/md01/nieyg/ref/hard-mask/genome_modify/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt" # same as before

nohup cellranger-arc mkref --config=hg38_hard_masked.config &

# 2. Cellranger (STAR) by masked genome

# rawdata pwd: /data/R03/zhangwx/project/human_PBMC/lib5_NhPBMC/
# ref: /md01/nieyg/ref/hard-mask/hg38_hard_masked

#PBS -N joint_FA
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome
cellranger-arc count --id=PBMC_4donor \
                       --reference=/md01/nieyg/ref/hard-mask/hg38_hard_masked \
                       --libraries=/data/R02/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/libraries.csv \
                       --localcores=24 \
                       --localmem=64
