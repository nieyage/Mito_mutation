
mut_data_high_con<- read.csv("mut_data_high_con_rescue_CA.csv")


# mutation signature

library(ggplot2)
library(dplyr)
library(stringr)

# 创建突变类型



# Step4: mutation signature 
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}
library(data.table)
# Process 3 digit signature based on letters
library(Biostrings)
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"
fasta_seqs <- readDNAStringSet(unshifted_chrM_ref)

# 查找chrM序列（根据实际序列名称调整）
seq_names <- names(fasta_seqs)
chrM_index <- grep("chrM|MT|mitochondria", seq_names, ignore.case = TRUE)
chrM_seq <- as.character(fasta_seqs[[chrM_index]])
ref_all <- data.frame(
  pos = 1:nchar(chrM_seq),
  ref = strsplit(chrM_seq, "")[[1]]
)
ref_all$ref <- toupper(ref_all$ref)
write.table(ref_all, "/md01/nieyg/ref/mito_ref/mm10/chrM_refAllele.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
random_cell_50<- sample(unique(mut_data_high_con$barcode),1)
subset<- mut_data_high_con[which(mut_data_high_con$barcode%in%random_cell_50),]
subset <- subset %>%
  mutate(
    mutation_type = paste0(Ref, ">", VarAllele)
  )
mutation<- paste0(subset$Position,subset$mutation_type)




called_variants <- mutation
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
  geom_bar(stat = "identity", position = "dodge") + #prettyGraphs::pretty_plot(fontsize = 8) + L_border() + 
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom",
    axis.text.x =element_text(angle=90,hjust=1,size=4)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y ="Substitution Rate\n(Expected / Observed)")+
  facet_wrap(~ group_change, scales = "free_x", nrow = 1)

pdf("./all_mito_signature_1cell.pdf", width = 8, height = 3)
p1;
dev.off()



