
nohup fastqc /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R1.fq.gz /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R2.fq.gz -o fastqc_results -t 16 &
nohup fastqc /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R1.fq.gz /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R2.fq.gz -o fastqc_results -t 16 &

#PBS -N seekarctools_unmasked
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/seekgene_multitissue/01_processed_data

seekarctools_py arc run \
    --rnafq1 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R1.fq.gz \
    --rnafq2 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R2.fq.gz \
    --atacfq1 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R1.fq.gz \
    --atacfq2 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R2.fq.gz \
    --samplename multitissue_test \
    --outdir /md01/nieyg/project/seekgene_multitissue/01_processed_data/seek_ref_unmasked \
    --refpath /md01/nieyg/ref/seekarc_ref/ensemble_102 \
    --include-introns \
    --core 24

# cellranger hard masked genome
# 当Cellranger-arc构建的参考基因组与SeekArcTools的STAR版本不兼容时，可以将cellranger-arc的STAR路径指定给SeekArcTools，例如：--star_path /path/to/cellranger-arc-2.0.2/lib/bin/STAR。

#PBS -N seekarctools_unmasked
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/seekgene_multitissue/01_processed_data

seekarctools_py arc run \
    --rnafq1 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R1.fq.gz \
    --rnafq2 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R2.fq.gz \
    --atacfq1 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R1.fq.gz \
    --atacfq2 /data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R2.fq.gz \
    --samplename multitissue_test \
    --outdir /md01/nieyg/project/seekgene_multitissue/01_processed_data/10x_ref_masked \
    --refpath /md01/nieyg/ref/hard-mask/mm10_hard_masked \
    --include-introns \
    --star_path /md01/nieyg/software/cellranger-arc-2.0.0/lib/bin/STAR
    --core 24







cd /md01/nieyg/project/seekgene_multitissue/01_processed_data
#zcat "/data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/rna/202511251801_A_SG251125-9_25050817_MT8_GE_L00_R2.fq.gz" | head -n $((1000 * 4)) > "subset.fastq"

seqtk sample -s "1234" "/data/R02/yuanwsy/multi_tissue/single_cell_MT8/seekgene/atac/202511251801_A_SG251125-9_25050817_MT8_arc_L00_R2.fq.gz" "1000" | gzip > "subset.fastq"
seqtk seq -A "subset.fastq" > sublset.fastq.fa
diamond blastx \
    --db /md01/nieyg/ncbi_database/nr.gz \
    --query sublset.fastq.fa \
    --out subset.fastq.fa.out \
    --outfmt 100 \
    --threads 1 \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --id 60 \
    --query-cover 70
diamond view \
    --daa subset.fastq.fa.out.daa \
    --out subset_matches.tsv

awk -F'\t' '{print $2}' "subset_matches.tsv" | \
    cut -d'|' -f2 | sort | uniq -c | sort -rn > "subset_matches_taxid_counts.txt"
echo "[$(date)] 步骤6: 获取物种名称..."
echo "TaxID,Count,ScientificName" > "subset_matches_species.csv"

./query_protein_info.sh



# libraries.csv
fastqs,sample,library_type
/data/R02/yuanwsy/multi_tissue/single_cell_MT8/10x/rna/20251129_LH00708_0261_B23FF52LT4/,MT8-RNA-251119,Gene Expression
/data/R02/yuanwsy/multi_tissue/single_cell_MT8/10x/atac/20251127_LH00708_0258_A23F5MLLT4/,MT8-ATAC-251119,Chromatin Accessibility


#PBS -N cellranger_masked
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=64000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/multitissue_plog/10x
cellranger-arc count --id=multitissue_10x_masked \
                       --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                       --libraries=libraries.csv \
                       --localcores=24 \
                       --localmem=64



# filter RNA and ATAC bam files 

nohup sh 1_filter_bam.sh &

# get fragments file from filtered bam file 
# 1. ATAC 
conda activate python38

newbam=atac_possorted_bam.dedup.bam
samtools index -@ 12 $newbam
sinto fragments -b $newbam -p 12 -f atac_fragments_filtered.tsv 
#sinto fragments -b filtered.bam -f fragments.tsv 

cat atac_fragments_filtered.tsv |sort -k1,1 -k2,2n | bgzip > atac_fragments_filtered.tsv.gz
tabix -b 2 -e 3 -p bed atac_fragments_filtered.tsv.gz

# 2. RNA 
# only keep the unique reads in bam file 
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/parse_cellsorted_bam.py .
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/run_dedup.sh .



nohup bash run_dedup.sh gex_possorted_bam &


    # 过滤配对好且MAPQ≥30的reads，解析CIGAR过滤匹配长度≥30
    samtools view -h -f 0x2 gex_possorted_bam.unique.bam | \
    awk 'BEGIN{OFS="\t"}
    /^@/ {print; next}
    {
        cigar=$6;
        m_sum=0;
        while(match(cigar, /([0-9]+)M/)) {
            m_sum += substr(cigar, RSTART, RLENGTH-1);
            cigar = substr(cigar, RSTART + RLENGTH);
        }
        if(m_sum >= 30) print $0;
    }' | \
    samtools view -b -o gex_possorted_bam_unique_filtered.bam -



