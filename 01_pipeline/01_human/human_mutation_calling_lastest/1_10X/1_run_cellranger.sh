# 2. Cellranger (STAR) by masked genome

# rawdata pwd: /data/R03/zhangwx/project/human_PBMC/lib5_NhPBMC/
# ref: /md01/nieyg/ref/hard-mask/hg38_hard_masked

#PBS -N 10_test
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/mito_mutation/01_pipeline/11_lastest
cellranger-arc count --id=PBMC_lib5 \
                       --reference=/md01/nieyg/ref/hard-mask/hg38_hard_masked_addflank \
                       --libraries=/data/R02/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/libraries.csv \
                       --localcores=32 \
                       --localmem=64
