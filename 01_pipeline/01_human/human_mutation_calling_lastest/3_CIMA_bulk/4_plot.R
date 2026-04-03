
# A. Coverage

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

coverage_file <- "./barcode_coverage.tsv.gz"
df <- read_tsv(coverage_file)

df <- df %>%
  mutate(
    position = as.integer(Position),
    coverage = as.numeric(Count)
  )

# 线粒体基因组特征区域
mt_regions <- data.frame(
  name = c("D-loop","D-loop", "12S rRNA", "16S rRNA", "ND1", "ND2", "CO1", "CO2", "ATP8", 
           "ATP6", "CO3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"),
  start = c(1,16024, 648, 1671, 3307, 4470, 5904, 7586, 8366, 8527, 9207, 
            10059, 10470, 10760, 12337, 14149, 14747),
  end = c(576,16569, 1601, 3229, 4262, 5511, 7445, 8269, 8572, 9207, 9990, 
          10404, 10766, 12137, 14148, 14673, 15887),
  type = c("D-loop","D-loop", "rRNA", "rRNA", "Protein", "Protein", "Protein", "Protein", 
           "Protein", "Protein", "Protein", "Protein", "Protein", "Protein", 
           "Protein", "Protein", "Protein")
)

plot_data <- df
max_coverage <- max(df$coverage, na.rm = TRUE)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)

p <- ggplot() +
  geom_line(data = plot_data, 
            aes(x = position, y = coverage, color = Barcode),
            size = 0.8) +
  scale_color_brewer(palette = "Set2") +
  geom_rect(data = mt_regions,
            aes(xmin = start, xmax = end, 
                ymin = max_coverage * 1.05, ymax = max_coverage * 1.1,
                fill = type),
            alpha = 0.7, inherit.aes = FALSE) +
  geom_text(data = mt_regions,
            aes(x = (start + end)/2, 
                y = max_coverage * 1.075, 
                label = name),
            size = 4, angle = 45, hjust = 1, vjust = 0.5, inherit.aes = FALSE) +
  labs(title = " ",
       x = "position (bp)", 
       y = "coverage",
       color = "celltype",
       fill = "Genomic Region") +
   scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 保存图形
output_file <- "mitochondrial_coverage_linearplot_percelltype.pdf"
ggsave(output_file, p, width = 14, height = 8)
message("图形已保存至: ", output_file)

# B. mtDNA % 
# 安装必要的包（如果尚未安装）
# BiocManager::install("Rsamtools")
# BiocManager::install("GenomicAlignments")

library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(data.table)

singlecell<- read.csv("/data/R01/renwx5/MT_mution/CIMA/align/new_ref/CNR1281138/outs/singlecell.csv")
celltype_file <- "/data/R01/renwx5/MT_mution/CIMA/annotation/CNR1281138_annotation.csv"
singlecell$mt_ratio<- (singlecell$mt_region_fragments/singlecell$fragments)*100
if (file.exists(celltype_file)) {
  celltype_df <- read_csv(celltype_file)
  
  # 检查必要的列是否存在
  if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
    celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
    message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")
    
    # 筛选有类型信息的细胞
    singlecell <- singlecell %>% filter(Cell %in% names(celltype_mapping))
    singlecell$celltype <- celltype_mapping[singlecell$Cell]
  } else {
    message("细胞类型文件缺少必要列，跳过细胞类型分析")
    singlecell$celltype <- "Unknown"
  }
} else {
  message("细胞类型文件不存在，跳过细胞类型分析")
  singlecell$celltype <- "Unknown"
}

library(ggplot2)
library(dplyr)

font_family <- "Arial"
p <- ggplot(singlecell, aes(x = celltype, y = mt_ratio, fill = celltype)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               color = "black", linewidth = 0.4, fatten = 1.5) +
  geom_point(position = position_jitter(width = 0.15, height = 0), 
             alpha = 0.4, size = 1.2, color = "black", shape = 16) +
  labs(
    x = "Cell Type",
    y = "mtDNA content (%)"
  ) +
  theme_classic() +
  theme(
    # 字体设置
    text = element_text(family = font_family, color = "black"),
    axis.text = element_text(family = font_family, color = "black", size = 10),
    axis.title = element_text(family = font_family, color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_brewer(palette = c("Set2"))
ggsave("./celltype_mtDNA_content_violin_enhanced.pdf", p, width = 11, height = 5.5, device = cairo_pdf)


# C. Mutation burden for each cell type
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# 步骤1: 读取变异稀疏矩阵数据
variant_file <- "./variant_sparse_matrix.tsv.gz"
message("正在读取变异数据: ", variant_file)

# 读取数据
variants <- fread(variant_file)
message("读取了 ", nrow(variants), " 行数据")
message("包含 ", length(unique(variants$cell)), " 个细胞")

# 查看数据结构
str(variants)
head(variants)

# 步骤2: 计算每个细胞中VAF不为0的mutation数目
message("计算每个细胞中VAF不为0的mutation数目...")

# 方法1: 直接计算
cell_mutations <- variants %>%
  # 筛选VAF > 0的突变（即alt_count > 0）
  filter(vaf > 0) %>%
  # 计算每个细胞的突变数
  group_by(cell) %>%
  summarise(
    mutation_count = n(),  # VAF > 0的突变数
    total_variants = n(),  # 同mutation_count
    mean_vaf = mean(vaf, na.rm = TRUE),
    max_vaf = max(vaf, na.rm = TRUE),
    .groups = "drop"
  )

# 按细胞类型统计
celltype_stats <- cell_mutations %>%
  group_by(cell) %>%
  summarise(
    n_cells = n(),
    mean_mutations = mean(mutation_count, na.rm = TRUE),
    median_mutations = median(mutation_count, na.rm = TRUE),
    sd_mutations = sd(mutation_count, na.rm = TRUE),
    max_mutations = max(mutation_count, na.rm = TRUE),
    cells_with_mutations = sum(mutation_count > 0),
    mutation_rate = cells_with_mutations / n_cells * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_mutations))

# 打印细胞类型统计
message("\n=== 细胞类型突变统计 ===")
print(celltype_stats)

# 计算每个细胞类型的均值和标准差
library(dplyr)
bar_data <- cell_mutations %>%
  group_by(cell) %>%
  summarise(mean_mut = mean(mutation_count, na.rm = TRUE),
            sd_mut = sd(mutation_count, na.rm = TRUE),
            .groups = "drop")

n_celltypes <- nrow(bar_data)

# 绘制条形图
p <- ggplot(bar_data, aes(x = cell, y = mean_mut, fill = cell)) +
  geom_col(width = 0.7, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_mut - sd_mut, ymax = mean_mut + sd_mut),
                width = 0.2, linewidth = 0.4, color = "black") +
  geom_text(aes(label = round(mean_mut, 2), y = mean_mut + sd_mut + 0.1),
            vjust = 0, size = 3.5, color = "black") +
  labs(title = "Mitochondrial Mutation Burden by Cell Type",
       x = "Cell Type",
       y = "Average Number of Mutations per Cell (VAF > 0)") +
  theme_classic() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 避免标签重叠
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_brewer(palette = "Set2")

# 保存图形
output_file <- "celltype_mutation_burden_barplot.pdf"
ggsave(output_file, p, width = 5, height = 6)
message("条形图已保存至: ", output_file)


# E. Mutation spectrum

library(ggplot2)
library(dplyr)
library(stringr)

somatic_snv<- "./final_mutation_report.tsv"
somatic_snv <- fread(somatic_snv)

mutation_list <- unique(paste0(somatic_snv$Position,paste(somatic_snv$Ref,somatic_snv$Alt,sep=">")))


# Step4: mutation signature 
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}
library(data.table)
# Process 3 digit signature based on letters
library(Biostrings)
unshifted_chrM_ref="/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta"
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
write.table(ref_all, "/md01/nieyg/ref/mito_ref/hg38/chrM_refAllele.txt", 
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
called_variants <- mutation_list
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

pdf("./all_somotic_mito_snv_signature_.pdf", width = 8, height = 3)
p1;
dev.off()

source("/md01/nieyg/project/mito_mutation/01_pipeline/10_v5/SNVcalling_test_v5/plot_mutation_signature.R")

# 使用示例：
# mutation_data <- data.frame(
#   Position = c(100, 200, 300, 100, 200, 400, 500, 100, 200),
#   Ref = c("A", "C", "G", "A", "C", "T", "A", "A", "C"),
#   VarAllele = c("G", "T", "A", "G", "T", "C", "G", "G", "T")
# )

mutation_data<- somatic_snv[,c("Position","Ref","Alt")]
colnames(mutation_data)<- c("Position","Ref","VarAllele")
plot_mut<- plot_and_save_mutation_signature(
   mutation_data,
   output_prefix = "PBMC_test",
   title_prefix = "PBMC_test",width = 8, height = 8
 )



# F. Mutation frequency spectrum
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

variant_file <- "variant_sparse_matrix.tsv.gz"
variants <- fread(variant_file)
variants <- variants[which(variants$vaf>0),]
variants$snv <- paste0(variants$pos,paste(variants$ref_base,variants$alt_base,sep=">"))
variants <- variants[which(variants$snv%in%mutation_list),]

variants_with_celltype <- variants %>%
  mutate(
    # 将VAF转换为百分比 (0-100)
    vaf_percent = vaf * 100,
    # 添加细胞类型
    #celltype = celltype_mapping[cell]
  ) 
# 步骤4: 创建频率分组
# 定义分组：0-10, 10-20, ..., 90-100
breaks <- seq(0, 100, by = 10)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1], "%")

variants_with_celltype <- variants_with_celltype %>%
  mutate(
    vaf_group = cut(vaf_percent, 
                    breaks = breaks, 
                    labels = labels,
                    include.lowest = TRUE, 
                    right = FALSE)
  )

# 查看分组情况
table(variants_with_celltype$vaf_group)


# 步骤5: 计算每个细胞类型在每个频率区间的突变数
freq_dist <- variants_with_celltype %>%
  filter(!is.na(vaf_group)) %>%  # 移除无法分组的记录
  group_by(cell, vaf_group) %>%
  summarise(
    mutation_count = n(),
    .groups = "drop"
  ) %>%
  # 确保所有细胞类型和所有分组都有值（填充0）
  complete(cell, vaf_group, fill = list(mutation_count = 0)) %>%
  # 按分组顺序排序
  mutate(vaf_group = factor(vaf_group, levels = labels))


# 分面图 - 每个细胞类型单独显示
p_facet <- ggplot(freq_dist, 
                  aes(x = vaf_group, y = mutation_count, fill = vaf_group)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black") +
  facet_wrap(~ cell, scales = "free_y", ncol = 3) +
  labs(
    title = "Mutation Frequency Distribution by Cell Type",
    x = "VAF Range (%)",
    y = "Number of Mutations"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

print(p_facet)
ggsave("mutation_freq_dist_facet.pdf", p_facet, width = 10, height = 10)
