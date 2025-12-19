
# A. Coverage

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# 步骤1: 读取覆盖度数据
message("步骤1: 读取覆盖度数据...")

# 从压缩文件读取数据
coverage_file <- "/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/barcode_coverage.tsv.gz"
df <- read_tsv(coverage_file, col_names = c("barcode", "position", "coverage"))
df <- df %>%
  mutate(
    position = as.integer(position),
    coverage = as.numeric(coverage)
  )

message("读取了 ", nrow(df), " 行数据")
message("包含 ", length(unique(df$barcode)), " 个细胞")

# 步骤2: 读取细胞类型信息
message("步骤2: 加载细胞类型信息...")

# 假设细胞类型信息文件路径
celltype_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"

if (file.exists(celltype_file)) {
  celltype_df <- read_csv(celltype_file)
  
  # 检查必要的列是否存在
  if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
    celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
    message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")
    
    # 筛选有类型信息的细胞
    df <- df %>% filter(barcode %in% names(celltype_mapping))
    df$celltype <- celltype_mapping[df$barcode]
  } else {
    message("细胞类型文件缺少必要列，跳过细胞类型分析")
    df$celltype <- "Unknown"
  }
} else {
  message("细胞类型文件不存在，跳过细胞类型分析")
  df$celltype <- "Unknown"
}

# 步骤3: 计算每个细胞类型的平均覆盖度
message("步骤3: 计算平均覆盖度...")

# 计算每个细胞类型的平均覆盖度
avg_coverage <- df %>%
  group_by(celltype, position) %>%
  summarise(avg_coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  complete(celltype, position = 1:16569, fill = list(avg_coverage = 0))

# 步骤4: 计算统计信息
message("步骤4: 计算统计信息...")

celltype_stats <- df %>%
  group_by(barcode, celltype) %>%
  summarise(
    total_coverage = sum(coverage),
    .groups = "drop"
  ) %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    mean_coverage = mean(total_coverage / 16569),
    .groups = "drop"
  )

# 打印统计信息
print(celltype_stats)

# 步骤5: 绘制覆盖度图
message("步骤5: 绘制覆盖度图...")

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
# 合并统计信息
plot_data <- avg_coverage %>%
  left_join(celltype_stats, by = "celltype") %>%
  mutate(celltype_label = paste0(celltype, "(n=", n_cells, ")"))

max_coverage <- max(plot_data$avg_coverage, na.rm = TRUE)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)

p <- ggplot() +
  geom_line(data = plot_data, 
            aes(x = position, y = avg_coverage, color = celltype_label),
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


# 使用samtools提取必要信息，然后用R处理
extract_bam_stats <- function(bam_file, barcode_file) {
  # 读取目标barcode
  target_barcodes <- fread(barcode_file, header = FALSE)$V1
  
  # 使用samtools提取信息
  message("使用samtools提取BAM信息...")
  cmd <- sprintf("samtools view %s | awk '{barcode=\"\"; for(i=12;i<=NF;i++) if($i~/^CB:Z:/) {barcode=substr($i,6); break} if(barcode!=\"\") print barcode\"\\t\"$3}'", bam_file)
  
  # 读取数据
  data <- fread(cmd = cmd, sep = "\t", header = FALSE, 
                col.names = c("barcode", "chr"))
  
  message(sprintf("提取了 %d 条记录", nrow(data)))
  
  # 过滤目标barcode
  filtered <- data[barcode %in% target_barcodes, ]
  message(sprintf("目标barcode的记录: %d 条", nrow(filtered)))
  
  # 统计
  total <- filtered[, .(total_reads = .N), by = barcode]
  mt <- filtered[chr %in% c("chrM", "MT", "M", "chrMT"), 
                .(mt_reads = .N), by = barcode]
  
  # 合并
  result <- merge(total, mt, by = "barcode", all.x = TRUE)
  result[is.na(mt_reads), mt_reads := 0]
  result[, mt_ratio := mt_reads / total_reads * 100]
  
  return(result[order(-total_reads)])
}

bam_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam"
output_file <- "barcode_mtDNA_stats_from_bam.csv"

# 运行
stats <- extract_bam_stats("/md01/nieyg/project/mito_mutation/01_pipeline/01_NhPBMC_joint_masked_genome/PBMC_4donor/outs/atac_possorted_bam.bam", "/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/barcodes.txt")
print(head(stats, 10))
fwrite(stats, "barcode_stats_final.csv")

stats<- read.csv("barcode_stats_final.csv")
celltype_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"

if (file.exists(celltype_file)) {
  celltype_df <- read_csv(celltype_file)
  
  # 检查必要的列是否存在
  if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
    celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
    message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")
    
    # 筛选有类型信息的细胞
    stats <- stats %>% filter(barcode %in% names(celltype_mapping))
    stats$celltype <- celltype_mapping[stats$barcode]
  } else {
    message("细胞类型文件缺少必要列，跳过细胞类型分析")
    stats$celltype <- "Unknown"
  }
} else {
  message("细胞类型文件不存在，跳过细胞类型分析")
  stats$celltype <- "Unknown"
}

library(ggplot2)
library(dplyr)

font_family <- "Arial"
p <- ggplot(stats, aes(x = celltype, y = mt_ratio, fill = celltype)) +
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
variant_file <- "/md01/jinxu/Project/mgatk-speedup/44_GenoByCell/variant_sparse_matrix.tsv"
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

# 方法2: 更详细的统计（包括总突变数）
cell_mutations_detailed <- variants %>%
  group_by(cell) %>%
  summarise(
    # 总变异位点数
    total_variant_sites = n(),
    # VAF > 0的突变数
    mutation_count = sum(vaf > 0, na.rm = TRUE),
    # VAF > 0.1的突变数（高频突变）
    high_vaf_mutations = sum(vaf > 0.1, na.rm = TRUE),
    # VAF > 0.5的突变数（优势突变）
    dominant_mutations = sum(vaf > 0.5, na.rm = TRUE),
    # 平均VAF（只计算VAF > 0的）
    mean_vaf = mean(vaf[vaf > 0], na.rm = TRUE),
    # 中位数VAF
    median_vaf = median(vaf[vaf > 0], na.rm = TRUE),
    .groups = "drop"
  )

# 查看结果
message("\n=== 突变统计摘要 ===")
cat("总细胞数:", nrow(cell_mutations), "\n")
cat("平均每个细胞的突变数:", mean(cell_mutations$mutation_count), "\n")
cat("有突变的细胞比例:", 
    sum(cell_mutations$mutation_count > 0) / nrow(cell_mutations) * 100, "%\n")
cat("最大突变数:", max(cell_mutations$mutation_count), "\n")

# 步骤3: 读取细胞类型信息
message("\n步骤3: 读取细胞类型信息...")

# 假设细胞类型信息文件路径
celltype_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"

if (file.exists(celltype_file)) {
  celltype_df <- fread(celltype_file)
  
  # 检查必要的列是否存在
  if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
    celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
    message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")
    
    # 合并细胞类型信息
    cell_mutations$celltype <- celltype_mapping[cell_mutations$cell]
    
    # 处理没有类型信息的细胞
    na_cells <- sum(is.na(cell_mutations$celltype))
    if (na_cells > 0) {
      message(na_cells, " 个细胞没有类型信息，标记为'Unknown'")
      cell_mutations$celltype[is.na(cell_mutations$celltype)] <- "Unknown"
    }
  } else {
    message("细胞类型文件缺少必要列，跳过细胞类型分析")
    cell_mutations$celltype <- "Unknown"
  }
} else {
  message("细胞类型文件不存在，跳过细胞类型分析")
  cell_mutations$celltype <- "Unknown"
}

# 步骤4: 数据清洗和准备
message("\n步骤4: 准备绘图数据...")

# 移除Unknown细胞类型（可选）
if ("Unknown" %in% unique(cell_mutations$celltype)) {
  original_cells <- nrow(cell_mutations)
  cell_mutations <- cell_mutations %>% filter(celltype != "Unknown")
  message("移除了 ", original_cells - nrow(cell_mutations), " 个Unknown细胞类型")
}

# 按细胞类型统计
celltype_stats <- cell_mutations %>%
  group_by(celltype) %>%
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

# 步骤5: 绘制小提琴图
message("\n步骤5: 绘制小提琴图...")

p <- ggplot(cell_mutations, aes(x = celltype, y = mutation_count, fill = celltype)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               color = "black", linewidth = 0.4, fatten = 1.5) +
  geom_point(position = position_jitter(width = 0.15, height = 0), 
             alpha = 0.4, size = 1.2, color = "black", shape = 16) +
  labs(title="Mitochondrial Mutation Burden by Cell Type",
    x = "Cell Type",
    y = "Number of Mutations per Cell (VAF > 0)"
  ) +
  theme_classic() +
  theme(
    # 字体设置
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text( color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_brewer(palette = c("Set2"))
# 保存图形
output_file <- "celltype_mutation_burden_violin.pdf"
ggsave(output_file, p, width = max(8, n_celltypes * 1.5), height = 6)
message("小提琴图已保存至: ", output_file)

# 7.2 突变数分布直方图
p_hist <- ggplot(cell_mutations, aes(x = mutation_count)) +
  geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = mean(cell_mutations$mutation_count), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Mutation Count per Cell",
    x = "Number of mutations (VAF > 0)",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave("mutation_count_distribution.pdf", p_hist, width = 8, height = 5)

# 步骤8: 保存统计结果
message("\n步骤8: 保存统计结果...")

# 保存详细的统计结果
output_stats_file <- "celltype_mutation_statistics.csv"
fwrite(celltype_stats, output_stats_file)
message("细胞类型统计已保存至: ", output_stats_file)

# 保存每个细胞的详细统计
output_cell_stats <- "cell_level_mutation_statistics.csv"
fwrite(cell_mutations, output_cell_stats)
message("细胞级别统计已保存至: ", output_cell_stats)

# D. Mutation burden normalized by mtDNA%

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# 检查数据
head(cell_mutations)
head(stats)

# 步骤1: 合并两个数据框
merged_data <- cell_mutations %>%
  rename(barcode = cell) %>%  # 统一列名
  inner_join(stats, by = "barcode") %>%
  # 确保没有重复
  distinct(barcode, .keep_all = TRUE)


# 步骤2: 计算normalized mutation burden
# 这里有几个可能的定义：
# 1. mutation_count / mt_ratio: 突变数除以mtDNA百分比
# 2. mutation_count / mt_reads: 突变数除以mtDNA reads数
# 3. mutation_count / (mt_reads/total_reads): 与1相同

merged_data <- merged_data %>%
  mutate(
    # 方法1: 直接用突变数除以mtDNA百分比
    mutation_burden_norm1 = ifelse(mt_ratio > 0, mutation_count / mt_ratio, NA),
    # 方法2: 用突变数除以mtDNA reads数（并乘以缩放因子，如1000）
    mutation_burden_norm2 = ifelse(mt_reads > 0, mutation_count / mt_reads * 1000, NA),
     )

# 步骤4: 按细胞类型统计
celltype_stats <- merged_data %>%
  group_by(celltype.x) %>%
  summarise(
    n_cells = n(),
    mean_mutation_count = mean(mutation_count, na.rm = TRUE),
    mean_mt_ratio = mean(mt_ratio, na.rm = TRUE),
    mean_norm_burden = mean(mutation_burden_norm1, na.rm = TRUE),
    median_norm_burden = median(mutation_burden_norm1, na.rm = TRUE),
    sd_norm_burden = sd(mutation_burden_norm1, na.rm = TRUE),
    correlation = cor(mutation_count, mt_ratio, use = "complete.obs", method = "spearman"),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_norm_burden))


# 步骤5: 可视化 - 按细胞类型分组的小提琴图
p <- ggplot(merged_data, aes(x = celltype.x, y = mutation_burden_norm1, fill = celltype.x)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               color = "black", linewidth = 0.4, fatten = 1.5) +
  geom_point(position = position_jitter(width = 0.15, height = 0), 
             alpha = 0.4, size = 1.2, color = "black", shape = 16) +
  labs(title="Mitochondrial Mutation Burden Normalized by mtDNA Content per Cell Type",
    x = "Cell Type",
    y = "Normlized Mutation Burden (mutations / mtDNA(%) )"
  ) +
  theme_classic() +
  theme(
    # 字体设置
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text( color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_brewer(palette = c("Set2"))
# 保存图形
output_file <- "mutation_burden_normalized_violin.pdf"
ggsave(output_file, p, width = max(8, n_celltypes * 1.5), height = 6)
message("小提琴图已保存至: ", output_file)


output_file <- "mutation_burden_normalized_data.csv"
fwrite(merged_data, output_file)
message("结果已保存至: ", output_file)

# 保存细胞类型统计
stats_file <- "celltype_normalized_burden_stats.csv"
fwrite(celltype_stats, stats_file)
message("细胞类型统计已保存至: ", stats_file)

# E. Mutation spectrum

library(ggplot2)
library(dplyr)
library(stringr)

somatic_snv<- "/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.somatic.tsv"
somatic_snv <- fread(somatic_snv)

colnames(somatic_snv)<- c("position","ref","alt","ref_fw","ref_rev","alt_fw","alt_rev","strand_score","mean_vaf","var_vaf","lis","pct_conf","pct_vaf_pos")
mutation_list <- unique(paste0(somatic_snv$position,paste(somatic_snv$ref,somatic_snv$alt,sep=">")))


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

pdf("./all_somotic_mito_snv_signature.pdf", width = 8, height = 3)
p1;
dev.off()

# F. Mutation frequency spectrum
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 步骤1: 读取原始变异数据
variant_file <- "/md01/jinxu/Project/mgatk-speedup/44_GenoByCell/variant_sparse_matrix.tsv"
variants <- fread(variant_file)
variants <- variants[which(variants$vaf>0),]
variants$snv <- paste0(variants$pos,paste(variants$ref_base,variants$alt_base,sep=">"))
variants <- variants[which(variants$snv%in%mutation_list),]

# 步骤2: 读取细胞类型信息
celltype_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"
celltype_df <- fread(celltype_file)
celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)

# 步骤3: 合并细胞类型信息到变异数据
variants_with_celltype <- variants %>%
  mutate(
    # 将VAF转换为百分比 (0-100)
    vaf_percent = vaf * 100,
    # 添加细胞类型
    celltype = celltype_mapping[cell]
  ) %>%
  # 移除没有细胞类型的记录
  filter(!is.na(celltype) & celltype != "")

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
  group_by(celltype, vaf_group) %>%
  summarise(
    mutation_count = n(),
    .groups = "drop"
  ) %>%
  # 确保所有细胞类型和所有分组都有值（填充0）
  complete(celltype, vaf_group, fill = list(mutation_count = 0)) %>%
  # 按分组顺序排序
  mutate(vaf_group = factor(vaf_group, levels = labels))


# 分面图 - 每个细胞类型单独显示
p_facet <- ggplot(freq_dist, 
                  aes(x = vaf_group, y = mutation_count, fill = vaf_group)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black") +
  facet_wrap(~ celltype, scales = "free_y", ncol = 3) +
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

# G. max VAF vs.  population score  

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

head(somatic_snv)

# 安装并加载ggpubr包
# install.packages("ggpubr")
library(ggpubr)

p1 <- ggplot(somatic_snv, aes(x = mean_vaf, y = pct_conf)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  # 使用stat_cor添加统计信息
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x.npc = "left", label.y.npc = "top",
           size = 5, color = "black", fontface = "bold") +
  labs(
    title = "Mean VAF vs Confidence Percentage",
    x = "Mean VAF (log10)",
    y = "Confidence Percentage (log10)",
    subtitle = paste("n =", nrow(somatic_snv), "mutations")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

p2 <- ggplot(somatic_snv, aes(x = mean_vaf, y = pct_vaf_pos)) +
  geom_point(alpha = 0.6, size = 2, color = "darkorange") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  # 使用stat_cor添加统计信息
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x.npc = "left", label.y.npc = "top",
           size = 5, color = "black", fontface = "bold") +
  labs(
    title = "Mean VAF vs VAF-positive Cells Percentage",
    x = "Mean VAF (log10)",
    y = "VAF-positive Cells % (log10)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# 保存
ggsave("mean_vaf_vs_pct_conf.pdf", p1, width = 8, height = 6)
ggsave("mean_vaf_vs_pct_vaf_pos.pdf", p2, width = 8, height = 6)

# H. Lineage informative score (LIS) distribution 

# 绘制 lis 值密度分布图
p_lis_density <- ggplot(somatic_snv, aes(x = lis)) +
  geom_density(fill = "steelblue", alpha = 0.6, color = "steelblue", linewidth = 1) +
  # 添加中位数和均值垂直线
  geom_vline(aes(xintercept = median(lis, na.rm = TRUE)), 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean(lis, na.rm = TRUE)), 
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  # 添加统计信息标签
  annotate("text",
           x = quantile(somatic_snv$lis, 0.98, na.rm = TRUE),  # 右侧位置
           y = max(density(somatic_snv$lis, na.rm = TRUE)$y) * 0.9,
           label = paste0("n = ", format(nrow(somatic_snv), big.mark = ","), "\n",
                          "Mean = ", round(mean(somatic_snv$lis, na.rm = TRUE), 5), "\n",
                          "Median = ", round(median(somatic_snv$lis, na.rm = TRUE), 5), "\n",
                          "SD = ", round(sd(somatic_snv$lis, na.rm = TRUE), 5), "\n",
                          "Min = ", round(min(somatic_snv$lis, na.rm = TRUE), 5), "\n",
                          "Max = ", round(max(somatic_snv$lis, na.rm = TRUE), 5)),
           hjust = 1, vjust = 1, size = 4,
           color = "black", fontface = "bold") +
  labs(
    title = "Distribution of LIS Values",
    x = "LIS (Likelihood Score)",
    y = "Density",
    subtitle = "All somatic SNV mutations"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.2)
  )

# 保存图形
ggsave("lis_density_distribution.pdf", p_lis_density, width = 10, height = 7)
print(p_lis_density)




