# A. coverage 

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(viridis)

count_directory <- arg_1
csv_file <- csv_file
# 步骤1: 加载细胞类型信息
message("步骤1: 加载细胞类型信息...")
celltype_df <- read.csv(csv_file)
celltype_mapping <- setNames(celltype_df$celltype, celltype_df$barcode)
message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")

# 步骤2: 解析count文件并计算每个细胞类型的覆盖度
message("步骤2: 解析count文件...")
count_files <- list.files(count_directory, pattern = "\\.(txt|counts)$", full.names = TRUE)
if (length(count_files) == 0) {
  count_files <- list.files(count_directory, full.names = TRUE)
}

message("找到 ", length(count_files), " 个count文件")

# 初始化存储结构
celltype_coverage <- list()
cell_stats <- data.frame()

for (file_path in count_files) {
  tryCatch({
    # 从文件名提取细胞ID
    cell_id <- tools::file_path_sans_ext(basename(file_path))
    
    if (!cell_id %in% names(celltype_mapping)) {
      next
    }
    
    celltype <- celltype_mapping[cell_id]
    
    # 读取count文件
    df <- read.table(file_path, sep = "\t", header = FALSE, 
                    col.names = c("chr", "pos", "base", "coverage", "A", "B", "C", "D", "E","F"))
    
    # 计算该细胞的统计信息
    total_coverage <- sum(df$coverage)
    zero_coverage_positions <- sum(df$coverage == 0)
    coverage_rate <- 1 - (zero_coverage_positions / nrow(df))
    
    cell_stats <- rbind(cell_stats, data.frame(
      cell_id = cell_id,
      celltype = celltype,
      total_coverage = total_coverage,
      zero_positions = zero_coverage_positions,
      coverage_rate = coverage_rate
    ))
    
    # 按位置排序
    df <- df[order(df$pos), ]
    
    # 确保覆盖所有位置 (1-16569)
    coverage_vector <- numeric(16569)
    for (i in 1:nrow(df)) {
      pos <- df$pos[i] - 1  # 转换为0-based索引
      if (pos >= 0 && pos < 16569) {
        coverage_vector[pos + 1] <- df$coverage[i]
      }
    }
    
    # 存储到对应细胞类型的列表
    if (is.null(celltype_coverage[[celltype]])) {
      celltype_coverage[[celltype]] <- list()
    }
    celltype_coverage[[celltype]] <- c(celltype_coverage[[celltype]], list(coverage_vector))
    
  }, error = function(e) {
    message("处理文件 ", file_path, " 时出错: ", e$message)
  })
}

# 步骤3: 计算每个细胞类型的平均覆盖度和污染率
message("步骤3: 计算细胞类型统计信息...")

# 计算每个细胞类型的平均覆盖度
avg_coverage <- list()
for (celltype in names(celltype_coverage)) {
  coverage_arrays <- celltype_coverage[[celltype]]
  if (length(coverage_arrays) > 0) {
    coverage_matrix <- do.call(rbind, coverage_arrays)
    avg_coverage[[celltype]] <- colMeans(coverage_matrix)
  } else {
    avg_coverage[[celltype]] <- numeric(16569)
  }
}

# 计算每个细胞类型的污染率（基于零覆盖位置的比例）
celltype_stats <- cell_stats %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    mean_total_coverage = mean(total_coverage),
    mean_coverage_rate = mean(coverage_rate),
    contamination_rate = 1 - mean_coverage_rate,  # 污染率 = 1 - 覆盖度比例
    sd_coverage_rate = sd(coverage_rate)
  ) %>%
  arrange(desc(mean_coverage_rate))

# 打印细胞类型统计信息
message("\n细胞类型统计信息:")
print(celltype_stats)

# 步骤4: 准备绘图数据
message("步骤4: 准备绘图数据...")

# 创建用于ggplot2的线性图数据
plot_data_linear <- data.frame()
for (celltype in names(avg_coverage)) {
  coverage <- avg_coverage[[celltype]]
  temp_df <- data.frame(
    position = 1:16569,
    coverage = coverage,
    celltype = celltype
  )
  plot_data_linear <- rbind(plot_data_linear, temp_df)
}

# 添加细胞类型统计信息到绘图数据
celltype_info <- celltype_stats %>%
  select(celltype, mean_coverage_rate, contamination_rate, n_cells)

plot_data_linear <- plot_data_linear %>%
  left_join(celltype_info, by = "celltype") %>%
  mutate(celltype_label = paste0(celltype, " (n=", n_cells))

# 步骤5: 使用ggplot2绘制线性视图
message("步骤5: 绘制线性视图...")
# 带颜色编码的版本
# 小鼠线粒体基因组特征区域（基于MM10/GRCm38）
mt_regions <- data.frame(
  name = c("D-loop", "tRNA-Phe", "12S rRNA", "tRNA-Val", "16S rRNA", "tRNA-Leu", 
           "ND1", "tRNA-Ile", "tRNA-Gln", "tRNA-Met", "ND2", "tRNA-Trp", 
           "tRNA-Ala", "tRNA-Asn", "tRNA-Cys", "tRNA-Tyr", "CO1", 
           "tRNA-Ser", "tRNA-Asp", "CO2", "tRNA-Lys", "ATP8", "ATP6", 
           "CO3", "tRNA-Gly", "ND3", "tRNA-Arg", "ND4L", "ND4", 
           "tRNA-His", "tRNA-Ser", "tRNA-Leu", "ND5", "ND6", "tRNA-Glu", 
           "CYTB", "tRNA-Thr", "tRNA-Pro", "D-loop"),
  start = c(15437, 1, 62, 957, 1025, 2012, 2071, 3055, 3118, 3192, 3260, 
            4209, 4278, 4345, 4413, 4480, 4544, 6109, 6176, 6243, 6936, 
            7007, 7070, 7749, 8438, 8505, 8889, 8965, 9060, 10365, 10434, 
            10501, 10570, 12040, 12706, 12772, 13856, 13923, 13990),
  end = c(16299, 61, 956, 1024, 2011, 2070, 3054, 3117, 3191, 3259, 4208, 
          4277, 4344, 4412, 4479, 4543, 6108, 6175, 6242, 6935, 7006, 
          7069, 7748, 8437, 8504, 8888, 8964, 9059, 10364, 10433, 10500, 
          10569, 12039, 12705, 12771, 13855, 13922, 13989, 15436),
  type = c("D-loop", "tRNA", "rRNA", "tRNA", "rRNA", "tRNA", "Protein", 
           "tRNA", "tRNA", "tRNA", "Protein", "tRNA", "tRNA", "tRNA", 
           "tRNA", "tRNA", "Protein", "tRNA", "tRNA", "Protein", "tRNA", 
           "Protein", "Protein", "Protein", "tRNA", "Protein", "tRNA", 
           "Protein", "Protein", "tRNA", "tRNA", "tRNA", "Protein", 
           "Protein", "tRNA", "Protein", "tRNA", "tRNA", "D-loop")
)

# 转换为数值类型
mt_regions$start <- as.numeric(mt_regions$start)
mt_regions$end <- as.numeric(mt_regions$end)

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")


# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)
max_coverage <- max(plot_data_linear$coverage, na.rm = TRUE)
p_linear_with_track <- ggplot() +
  geom_line(data = plot_data_linear, 
            aes(x = position, y = coverage, color = celltype_label),
            size = 0.8) +
  scale_color_manual(values = myUmapcolors) +
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

linear_output <- "./mitochondrial_coverage_linear_with_contamination.pdf"
ggsave(linear_output, p_linear_with_track, width = 16, height = 8)
message("线性图已保存至: ", linear_output)






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

bam_file <- "/md01/nieyg/project/mito_mutation/02_mm10_pipeline/01_masked_cellranger/BMMC_27m_ATAC/outs/atac_possorted_bam.bam"
barcode_file<- "/md01/nieyg/project/mito_mutation/02_mm10_pipeline/02_modify_bam/barcodes.txt"
output_file <- "barcode_mtDNA_stats_from_bam.csv"

# 运行
stats <- extract_bam_stats(bam_file, barcode_file)
print(head(stats, 10))
fwrite(stats, output_file)
csv_file="/md01/nieyg/project/mito_mutation/02_mm10_pipeline/03_singlecell_SNV/gexbarcode_celltype.csv" 
stats<- read.csv("barcode_stats_final.csv")
celltype_file <- csv_file


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


# C. mutation burden 

library(ggplot2)
library(dplyr)

# 计算每个细胞的突变数目
cell_mutation_counts <- mut_data_high_con %>%
  group_by(barcode, celltype) %>%
  summarise(mutation_count = n(), .groups = "drop")

# 查看统计摘要
print("每个细胞类型的突变数目统计:")
print(cell_mutation_counts %>%
  group_by(celltype) %>%
  summarise(
    mean_count = mean(mutation_count),
    median_count = median(mutation_count),
    min_count = min(mutation_count),
    max_count = max(mutation_count),
    n_cells = n()
  ))

# 专业美化版本
font_family <- "Arial"
p_enhanced <- ggplot(cell_mutation_counts, aes(x = celltype, y = mutation_count, fill = celltype)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               color = "black", linewidth = 0.4, fatten = 1.5) +
  geom_point(position = position_jitter(width = 0.15, height = 0), 
             alpha = 0.4, size = 0.1, color = "black", shape = 16) +
  scale_color_manual(values = myUmapcolors) +
  scale_fill_manual(values = myUmapcolors)+
  labs(
    x = "Cell Type",
    y = "Number of Mutations per Cell"
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
  ) 
# 保存增强版
ggsave("./celltype_mutation_violin_enhanced.pdf", p_enhanced, width = 16, height = 5.5, device = cairo_pdf)

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


# E. mutation signature

library(ggplot2)
library(dplyr)
library(stringr)

# 创建突变类型
SNV_filter <- SNV_filter %>%
  mutate(
    mutation_type = paste0(Ref, ">", VarAllele)
  )
mutation<- paste0(SNV_filter$Position,SNV_filter$mutation_type)
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

pdf("./all_mito_signature.pdf", width = 8, height = 3)
p1;
dev.off()


# F.  Mutation frequency spectrum
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

variants_with_celltype <- mut_data_high_con %>%
  mutate(
    # 将VAF转换为百分比 (0-100)
    vaf_percent = VarFreq_numeric * 100
  ) %>%
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
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
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
ggsave("mutation_freq_dist_facet.pdf", p_facet, width = 15, height = 10)






