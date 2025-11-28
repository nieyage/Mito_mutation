# Step1: split bam by sample 
# Sample ID from Mitosort:/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv

#!/bin/bash

# 输入文件
csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"  # 请替换为实际的CSV文件路径
unshifted_bam="/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_unshifted.realign.bam"
shifted_bam="/md01/nieyg/project/mito_mutation/01_pipeline/02_getchrM_unmapped/PBMC_chrM_unmapped_bwa_shifted.realign.bam"
subset_bam_tool="/md01/nieyg/software/subset-bam_linux"

# 输出目录
output_dir="./donor_bams/"
mkdir -p "$output_dir"

# 临时目录存放barcode列表
barcode_dir="./barcode_lists/"
mkdir -p "$barcode_dir"

# 从CSV文件中提取每个Donor的barcode
echo "Extracting barcodes for each donor from CSV file..."

# 提取Donor1的barcode
awk -F',' '$2 == "\"Donor1\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor1_barcodes.txt"
echo "Donor1: $(wc -l < ${barcode_dir}donor1_barcodes.txt) barcodes"

# 提取Donor2的barcode
awk -F',' '$2 == "\"Donor2\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor2_barcodes.txt"
echo "Donor2: $(wc -l < ${barcode_dir}donor2_barcodes.txt) barcodes"

# 提取Donor3的barcode
awk -F',' '$2 == "\"Donor3\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor3_barcodes.txt"
echo "Donor3: $(wc -l < ${barcode_dir}donor3_barcodes.txt) barcodes"

# 提取Donor4的barcode
awk -F',' '$2 == "\"Donor4\"" {gsub(/"/, "", $1); print $1}' "$csv_file" > "${barcode_dir}donor4_barcodes.txt"
echo "Donor4: $(wc -l < ${barcode_dir}donor4_barcodes.txt) barcodes"

# 为每个Donor提取BAM文件
echo "Extracting BAM files for each donor..."

process_donor() {
    local donor=$1
    local barcode_file="${barcode_dir}${donor}_barcodes.txt"
    
    echo "Processing $donor..."
    
    # 提取unshifted BAM
    echo "  Extracting unshifted BAM for $donor..."
    "$subset_bam_tool" --bam "$unshifted_bam" --cell-barcodes "$barcode_file" --cores 1 --out-bam "${output_dir}${donor}_unshifted.bam"
    
    # 提取shifted BAM
    echo "  Extracting shifted BAM for $donor..."
    "$subset_bam_tool" --bam "$shifted_bam" --cell-barcodes "$barcode_file" --cores 1 --out-bam "${output_dir}${donor}_shifted.bam"
    
    echo "  Finished processing $donor"
}

# 处理所有Donor
process_donor "donor1"
process_donor "donor2"
process_donor "donor3"
process_donor "donor4"

echo "All Done!"
echo "Output files are in: $output_dir"


# Step2: Call germline mutation for each donor 

#!/bin/bash
# 设置变量
unshifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.fasta
shifted_chrM_ref=/md01/nieyg/ref/mito_ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
chrM_len="/md01/nieyg/ref/hard-mask/hg38_hard_masked/hg38_chrM.len"
picard_jar="/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
pileup2snp_tool="varscan pileup2snp"  # 请确保这个工具在PATH中
pileup_inf_tool="/md01/jinxu/bin/pileup_inf_rj.pl"

# donor 列表
donors=("donor1" "donor2" "donor3" "donor4")

# 处理每个 donor
for donor in "${donors[@]}"; do
    echo "Processing $donor..."
    
    # 创建 donor 目录（如果不存在）
    mkdir -p "$donor"
    
    # 定义文件路径
    unshifted_bam="./donor_bams/${donor}_unshifted.bam"
    shifted_bam="./donor_bams/${donor}_shifted.bam"
    
    # 检查输入文件是否存在
    if [[ ! -f "$unshifted_bam" ]]; then
        echo "Error: Unshifted BAM file not found: $unshifted_bam"
        continue
    fi
    if [[ ! -f "$shifted_bam" ]]; then
        echo "Error: Shifted BAM file not found: $shifted_bam"
        continue
    fi
    
    # 1. 处理 unshifted BAM 文件
    echo "  Processing unshifted BAM..."
    
    # 排序
    samtools sort "$unshifted_bam" -o "${donor}/${donor}_unshifted.sorted.bam"
    
    # 标记重复
    java -Xmx4g -jar "$picard_jar" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${donor}/${donor}_unshifted.sorted.bam" \
        OUTPUT="${donor}/${donor}_unshifted.rmdup.bam" \
        METRICS_FILE="${donor}/${donor}_unshifted.metrics"
    
    # 生成 mpileup
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$unshifted_chrM_ref" \
        -x "${donor}/${donor}_unshifted.rmdup.bam" > "${donor}/${donor}_unshifted.rmdup.mpileup"
    
    # 2. 处理 shifted BAM 文件
    echo "  Processing shifted BAM..."
    
    # 排序
    samtools sort "$shifted_bam" -o "${donor}/${donor}_shifted.sorted.bam"
    
    # 标记重复
    java -Xmx4g -jar "$picard_jar" \
        MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
        INPUT="${donor}/${donor}_shifted.sorted.bam" \
        OUTPUT="${donor}/${donor}_shifted.rmdup.bam" \
        METRICS_FILE="${donor}/${donor}_shifted.metrics"
    
    # 生成 mpileup
    samtools mpileup \
        -l "$chrM_len" \
        -q 30 -Q 30 -f "$shifted_chrM_ref" \
        -x "${donor}/${donor}_shifted.rmdup.bam" > "${donor}/${donor}_shifted.rmdup.mpileup"
    
    # 3. 合并 mpileup 文件
    echo "  Merging mpileup files..."
    
    shifted_mpileup="${donor}/${donor}_shifted.rmdup.mpileup"
    unshifted_mpileup="${donor}/${donor}_unshifted.rmdup.mpileup"
    output_mpileup="${donor}/${donor}_combined.mpileup"
    
    # 检查输入文件是否存在
    if [[ ! -f "$shifted_mpileup" ]]; then
        echo "Error: shifted mpileup file not found: $shifted_mpileup"
        continue
    fi
    if [[ ! -f "$unshifted_mpileup" ]]; then
        echo "Error: unshifted mpileup file not found: $unshifted_mpileup"
        continue
    fi
    
    # 临时文件
    temp_part1="${donor}/temp_part1.mpileup"
    temp_part2="${donor}/temp_part2.mpileup"
    temp_part3="${donor}/temp_part3.mpileup"
    
    # 清理函数
    cleanup() {
        rm -f "$temp_part1" "$temp_part2" "$temp_part3"
    }
    trap cleanup EXIT
    
    # 第一部分：从shifted文件中提取8570-9144行，位置改为1-575
    awk -v start=8570 -v end=9144 -v new_start=1 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        new_pos = $2 - start + new_start;
        $2 = new_pos;
        print $0;
    }' "$shifted_mpileup" > "$temp_part1"
    
    # 第二部分：从unshifted文件中提取576-16024行，位置不变
    awk -v start=576 -v end=16024 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        print $0;
    }' "$unshifted_mpileup" > "$temp_part2"
    
    # 第三部分：从shifted文件中提取8025-8569行，位置改为16025-16569
    awk -v start=8025 -v end=8569 -v new_start=16025 '
    BEGIN{OFS="\t"} $2 >= start && $2 <= end {
        new_pos = $2 - start + new_start;
        $2 = new_pos;
        print $0;
    }' "$shifted_mpileup" > "$temp_part3"
    
    # 合并三个部分
    cat "$temp_part1" "$temp_part2" "$temp_part3" > "$output_mpileup"
    
    # 清理临时文件
    # cleanup
    trap - EXIT
    
    # 验证合并结果
    total_lines=$(wc -l < "$output_mpileup")
    echo "  Combined mpileup lines: $total_lines"
    
    # 4. SNV calling 和计数
    echo "  Performing SNV calling and counting..."
    
    # SNV calling
    $pileup2snp_tool "$output_mpileup" --min-var-freq 0.01 --min-reads2 3 > "${donor}/${donor}.snv"
    
    # 计数
    $pileup_inf_tool "$output_mpileup" > "${donor}/${donor}.counts"
    
    echo "  Finished processing $donor"
    echo "----------------------------------------"
done

echo "All donors processed successfully!"


# for donors

library(dplyr)
donor1_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor1/donor1.snv", header = TRUE)
germline <- filter(donor1_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor1_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor2_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor2/donor2.snv", header = TRUE)
germline <- filter(donor2_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor2_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor3_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor3/donor3.snv", header = TRUE)
germline <- filter(donor3_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor3_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)

donor4_snv = read.table(file = "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/unmasked/donor4/donor4.snv", header = TRUE)
germline <- filter(donor4_snv, Reads2/(Reads1 + Reads2) > 0.9)
write.table(germline, file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor4_q20Q30.germline.snv",quote = FALSE,sep = "\t",row.names = F)


# old pipeline

library(data.table)

arg_1 = "/md01/nieyg/project/mito_mutation/01_pipeline/old_03_split_bam/unmasked_SNVcalling_percell/"

snv_files <- list.files(arg_1, pattern = "\\.snv$", full.names = TRUE)

# 读取并合并所有文件
mut_data <- rbindlist(lapply(snv_files, function(file) {
  # 读取数据
  dt <- fread(file)
  
  # 从文件名提取barcode（根据您的文件名格式调整）
  barcode <- gsub("\\.snv$", "", basename(file))
  
  # 创建突变标识符，添加barcode列
  dt[, `:=`(
    mutation_id = paste(Position, Ref, VarAllele, sep = "_"),
    barcode = barcode
  )]
  
  # 转换VarFreq为数值（去掉百分号）
  dt[, VarFreq_numeric := as.numeric(gsub("%", "", VarFreq)) / 100]
  
  return(dt)
}))

# 转换为宽格式矩阵
mut_matrix <- dcast(mut_data, 
                   barcode ~ mutation_id, 
                   value.var = "VarFreq_numeric",
                   fill = 0)

# 设置行名
final_matrix <- as.matrix(mut_matrix[, -1])
rownames(final_matrix) <- mut_matrix$barcode


germline1 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor1_q20Q30.germline.snv", header = T)
germline2 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor2_q20Q30.germline.snv", header = T)
germline3 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor3_q20Q30.germline.snv", header = T)
germline4 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor4_q20Q30.germline.snv", header = T)
germline1 <- paste(germline1$Position, germline1$Ref, germline1$VarAllele, sep = "_")
germline2 <- paste(germline2$Position, germline2$Ref, germline2$VarAllele, sep = "_")
germline3 <- paste(germline3$Position, germline3$Ref, germline3$VarAllele, sep = "_")
germline4 <- paste(germline4$Position, germline4$Ref, germline4$VarAllele, sep = "_")

germline_mutation<- c(germline1,germline2,germline3,germline4)



# 2.sc_SNV filter
# define file path and an empty dataframe
library(dplyr)
arg_1 = "/md01/nieyg/project/mito_mutation/01_pipeline/old_03_split_bam/unmasked_SNVcalling_percell/"
files <- dir(arg_1, pattern = "snv$")
path <- arg_1
i <- 1
final <- data.frame("Chrom" = 0 , "Position" = 0, "Ref" = 0, "VarAllele" = 0)
final <- final[-length(final$Chrom),]
sc_germline <- final


for(i in 1 : length(files)) {
  x <- paste0(path, files[i])
  y <- read.table(file = x, header = T, colClasses = c("character"))
  y$Reads1 <- as.numeric(y$Reads1)
  y$Reads2 <- as.numeric(y$Reads2)
  y$Reads2Plus <- as.numeric(y$Reads2Plus)
  y$Reads2Minus <- as.numeric(y$Reads2Minus)
  z <- filter(y, Reads2Plus > 1 & 
                Reads2Minus > 1 &
                Reads2Plus/(Reads2Minus + Reads2Plus) < 0.7 &
                Reads2Plus/(Reads2Minus + Reads2Plus) > 0.3 &
                (Reads1 + Reads2) >= 10 &
                Reads2/(Reads1 + Reads2) >= 0.10 &
                !(y$Ref == "G" & y$VarAllele == "T") &
                !(y$Ref == "C" & y$VarAllele == "A")
  )
  k <- filter(z, Reads2/(Reads1 + Reads2) >= 0.90)
  n <- select(z, c("Chrom", "Position", "Ref", "VarAllele"))
  k <- select(k, c("Chrom", "Position", "Ref", "VarAllele"))
  final <- merge(final, n, all = T)
  # final <- rbind(final, n)
  sc_germline <- rbind(sc_germline, k)
  # sc_germline <- merge(sc_germline, k, all = T)
}

frequency <- as.data.frame(table(sc_germline$Position))
frequency <- filter(frequency, frequency$Freq > length(files)*0.9)
final <- unique(final)
final_remove <- c()
for (i in 1:length(final$Position)) {
  if (final$Position[i] %in% frequency$Var1){
    final_remove <- append(final_remove, i)
  }
}
if(length(final_remove >= 1)){
  final <- final[-final_remove,]
}

write.table(final, file = "./MPP_cell_realign/MPP_cell_realign.sc.filter.snv", sep = "\t", quote = F, row.names = F)

# input germline and blacklist

germline1 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor1_q20Q30.germline.snv", header = T)
germline2 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor2_q20Q30.germline.snv", header = T)
germline3 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor3_q20Q30.germline.snv", header = T)
germline4 <- read.table(file = "/data/R02/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/donor4_q20Q30.germline.snv", header = T)
germline<- rbind(germline1,germline2,germline3,germline4)

blacklist <- read.table(file = "./blacklist.txt", header = T)

# remove germline mutation
SNV_remove <- c()
SNV_filter<- final
for(i in 1:length(SNV_filter$Position)) {
  if(SNV_filter$Position[i] %in% germline$Position) {
    SNV_remove <- append(SNV_remove, i)
  }
}
if(length(SNV_remove) >= 1){
  SNV_filter <- SNV_filter[-SNV_remove,]}

# remove mutation in black list
blacklist_pos <-c("302","309","311","312","313","316","514","515","523","524","3106","3109","3110")
blacklist_remove<- c()
for(i in 1:length(SNV_filter$Position)) {
  if(SNV_filter$Position[i] %in% blacklist_pos) {
    blacklist_remove <- append(blacklist_remove, i)
  }
}
if(length(blacklist_remove) >= 1){
  SNV_filter <- SNV_filter[-blacklist_remove,]}

# arrange
SNV_filter <- arrange(SNV_filter, Position)

high_con_mutation<- paste(SNV_filter$Position, SNV_filter$Ref, SNV_filter$VarAllele, sep = "_")


csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"  # 请替换为实际的CSV文件路径
metadata<- read.csv(csv_file)
mut_data$celltype<- metadata[match(mut_data$barcode,metadata$barcode),]$Annotation
mut_data$donor<- metadata[match(mut_data$barcode,metadata$barcode),]$sample

mut_data_high_con<- mut_data[which(mut_data$mutation_id%in%high_con_mutation),]


mut_data_high_con[order(mut_data_high_con$barcode),]


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
p_enhanced <- ggplot(cell_mutation_counts, aes(x = celltype, y = mutation_count, fill = celltype)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               color = "black", linewidth = 0.4, fatten = 1.5) +
  geom_point(position = position_jitter(width = 0.15, height = 0), 
             alpha = 0.4, size = 1.2, color = "black", shape = 16) +
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
  ) +
  scale_fill_brewer(palette = c("Set2"))
# 保存增强版
ggsave("../celltype_mutation_violin_enhanced.pdf", p_enhanced, width = 11, height = 5.5, device = cairo_pdf)


# Coverage 
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(viridis)

count_directory <- "/md01/nieyg/project/mito_mutation/01_pipeline/old_03_split_bam/unmasked_SNVcalling_percell"
csv_file <- "/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"

# 步骤1: 加载细胞类型信息
message("步骤1: 加载细胞类型信息...")
celltype_df <- read.csv(csv_file)
celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
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

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)
max_coverage <- max(plot_data_linear$coverage, na.rm = TRUE)
p_linear_with_track <- ggplot() +
  geom_line(data = plot_data_linear, 
            aes(x = position, y = coverage, color = celltype_label),
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

linear_output <- "../mitochondrial_coverage_linear_with_contamination.pdf"
ggsave(linear_output, p_linear_with_track, width = 14, height = 8)
message("线性图已保存至: ", linear_output)

# mutation signature

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

pdf("../all_mito_signature.pdf", width = 8, height = 3)
p1;
dev.off()

cowplot::ggsave2(p1, file = )

# for one barcode 




