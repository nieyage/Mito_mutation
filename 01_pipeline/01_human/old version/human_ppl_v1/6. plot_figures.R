# 1. The number of somatic mtDNA mutations per cell. (group by celltype)

library(dplyr)
final<-read.table(".../Data for plotting/cell_mutation_count.txt",row.names=1,header = T,sep="\t")
final$celltypes=factor(final$celltypes,levels=(c("HSC","LMPP","CLP","proB","preB","B","T","NK","CMP","GMP","MDP","pDC","mono","Ery")))

p <- ggplot(final, aes(x=celltypes,y=mucount)) +     
  geom_violin(aes(fill=celltypes),width = 1.1)
p <- p + theme_bw() + scale_y_continuous(limits =c(0,15),breaks = seq(0,15,2)) + 
  theme(
    panel.grid = element_blank(),  
    panel.border = element_blank(), 
    axis.line = element_line(size =1.2,colour = "black"), # 
    axis.ticks.length=unit(1.2, "pt"),  # 
    axis.ticks = element_line(size =1.2,colour = "black"), 
    axis.text.y = element_text(face="bold",size = 12,colour="black"),
    axis.text.x = element_text(face="bold",size = 12,colour="black"),
    axis.title.x = element_text(face="bold",size = 12),
    axis.title.y = element_text(size =12,face="bold"),
    legend.title = element_text(face="bold",size=12),
    legend.text=element_text(face="bold",size=12))+
labs(x="",y="mutation counts",fill="Group",title="mutation counts for all cell tyeps") +  
  guides(fill=F)  
p
ggsave(p ,filename = ".../mutation_percells.pdf",width = 20,height = 5)


# %mt reads.

p <- ggplot(plot_data, aes(x=stage,y=g_num)) +     
  geom_boxplot(aes(fill=stage))  
p <- p + theme_bw() + 
  scale_y_continuous(limits =c(0,0.8),breaks = seq(0,0.8,0.1)) +
  theme(
    panel.grid = element_blank(),   
    panel.border = element_blank(),  
    axis.line = element_line(size =1.2,colour = "black"), 
    axis.ticks.length=unit(2, "pt"),  
    axis.ticks = element_line(size =1.2,colour = "black"),  
    axis.text.y = element_text(face="bold",size = 12,colour="black"), 
    axis.text.x = element_text(face="bold",size = 12,colour="black"),
    axis.title.x = element_text(face="bold",size = 12),
    axis.title.y = element_text(size =12,face="bold"),
    legend.title = element_text(face="bold",size=12),
    legend.text=element_text(face="bold",size=12))+
labs(x="",y="%mt reads",fill="Group") +  
  guides(fill=F)  
p
ggsave(p,filename = ".../mt_percentage_mergecell.pdf",width = 20,height = 5)

# 3. Allele frequency spectrum of somatic mtDNA mutations. 


# plot the mutation signature 
# Any neighbouring SNVs will be merged into DBS/MBS variants.
library(MutationalPatterns)
library(BSgenome)
head(available.genomes())
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)
vcf_files <- list.files(path="/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/common_mutation/vcf/",pattern = ".vcf", full.names = TRUE)
sample_names <- c(  "AR3-C4", "AR3-C5")

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,group="circular", type = c("snv"))
time <- c( "AR3-C4", "AR3-C5")
#snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)
library("gridExtra")
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)
pdf("./AR3_pseudobulk_96_profile.pdf",width=12,height=6)
plot_96_profile(mut_mat)
dev.off()

# plot the location signature 
# 1. coverage 
AR3_C4_data<- read.table("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out/AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count")
AR3_C5_data<- read.table("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mito_sort_out/AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count")
AR3_C4_data$sample<-"AR3_C4"
AR3_C5_data$sample<-"AR3_C5"

mdf<- rbind(AR3_C4_data,AR3_C5_data)
library(ggplot2)
pdf("./AR3_pseudobulk_location_coverage_profile.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + 
  theme(legend.position = "none")

dev.off()

mdf$mutation<- (mdf$V9+mdf$V10)
pdf("./AR3_pseudobulk_location_alt_profile.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + 
  theme(legend.position = "none")

ggplot(mdf, aes(x = V2, y = Freq,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + 
  theme(legend.position = "none")


dev.off()


pdf("./AR3_pseudobulk_location_alt_profile16000.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + xlim(16080,16150)+
  theme(legend.position = "none")

ggplot(mdf, aes(x = V2, y = Freq)) + 
  geom_bar(stat="identity") +  
  #pretty_plot(fontsize = 8)  + 
  #scale_color_fill(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + xlim(16080,16150)+
  theme(legend.position = "none")
dev.off()




# Step1: QC metric for mutation 

Idents(combined)<- combined$detail_anno
combined$mito_reads_rate<- (combined$mitochondrial/combined$total)*100
library(scCustomize)
pdf("./09_mgatk_for_clean_bam/all_cell_type_mtDNA_depth.pdf",width=20,height=10)
VlnPlot(combined, c("mito_reads_rate","mtDNA_depth","nCount_alleles","nFeature_alleles"),col=myUmapcolors, pt.size = 0,ncol=1) + scale_y_log10()
VlnPlot(combined,log = TRUE, c("mito_reads_rate","mtDNA_depth","nCount_alleles","nFeature_alleles"), pt.size = 0,ncol=1) + scale_y_log10()
Stacked_VlnPlot(seurat_object = combined, features = c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), x_lab_rotate = TRUE,colors_use = myUmapcolors)
dev.off()

# Step2: VAF distribution
# global VAF distribution of mutation 
# the mutation freq barplot 
DefaultAssay(combined)<- "alleles"
data<- as.numeric(GetAssayData(combined))
# all mutation freq in single cell not all cells(remove the 0)
data<- data[data!=0]
data<- data.frame(variance=data)
data$variance<- data$variance;
pdf("./09_mgatk_for_clean_bam/All_celltype_mutation_freq.pdf",width=5,height=3)
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.1,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.1))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title="All_celltype_mutation_freq")+
  scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2))
  p
dev.off()

# VAF distribution of mutation in difference celltype 
# the mutation freq barplot 
Idents(combined)<- combined$Annotation
pdf("./09_mgatk_for_clean_bam/celltype_mutation_freq.pdf",width=10,height=6)
for (j in levels(combined)){
  print(j);
  obj<- subset(combined,idents=j)
  variable.sites <- IdentifyVariants(obj, assay = "mito", refallele = refAllele)
    # select the cutoff 
  cutoff_n_cells_conf_detected<- as.matrix(summary(variable.sites$n_cells_conf_detected))[2]
  cutoff_mean_coverage<- as.matrix(summary(variable.sites$mean_coverage))[2]
  if(cutoff_mean_coverage>40){cutoff_mean_coverage=40}else{cutoff_mean_coverage=cutoff_mean_coverage}
  if(cutoff_n_cells_conf_detected<5){cutoff_n_cells_conf_detected=5}else{cutoff_n_cells_conf_detected=cutoff_n_cells_conf_detected}
  
  # step1: Consider strand balance(DNA damage casued C-->A and G-->T error) and the VMR(rm the germline mutation) 
  # Step2: Restriction of minimum frequency and sequening depth:
  high.conf <- subset(
    variable.sites, 
    subset = n_cells_conf_detected >= cutoff_n_cells_conf_detected &
      strand_correlation >= 0.65 &
      mean_coverage >= cutoff_mean_coverage &
      vmr > 0.01
  )
  p1 <- VariantPlot(variants = variable.sites,concordance.threshold=0.65)
  print(p1)
  if(nrow(high.conf)>1){
  crc <- AlleleFreq(
    object = obj,
    variants = high.conf$variant,
    assay = "mito"
  )
  DefaultAssay(crc) <- "alleles"
  data<- GetAssayData(crc)
  d_varition1=as.vector(data)
  d_varition1=d_varition1[which(d_varition1!=0)]
  data<- data.frame(variance=d_varition1)
  data$variance<- data$variance;
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.1,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.1))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title=j)+
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=18,colour="black", face = "bold"),          
        axis.text.y = element_text(size=18,colour="black", face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"), 
        axis.title.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        line = element_line(size=1),
        axis.line = element_line(size =1.0,colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2))
    print(p)
}}
dev.off()

# Step3: mutation lacation 

# plot the location signature 
library(ggplot2)

high.conf<- read.csv("./09_mgatk_for_clean_bam/high_conf_alleles_info.csv")
head(high.conf)
tmp_data<- as.data.frame(table(high.conf$position))

pdf("./09_mgatk_for_clean_bam/AR3_pseudobulk_location_alt_profile.pdf",width=24,height=6)

ggplot(tmp_data, aes(x = Var1, y = Freq)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + 
  theme(legend.position = "none")
dev.off()

# Step4: mutation signature 
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}
library(data.table)
# Process 3 digit signature based on letters
ref_all <- fread("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/final/chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants
called_variants <- high.conf$variant
ref_all_long$called <- ref_all_long$variant %in% called_variants

# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
library(dplyr)
library(tidyr)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

# Visualize
library(tidyverse)
library(prettyGraphs)
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom",
    axis.text.x =element_text(angle=90,hjust=1,size=2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "./09_mgatk_for_clean_bam/all_mito_signature.pdf", width = 4, height = 2.4)

