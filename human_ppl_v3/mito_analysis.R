#!/usr/bin/env Rscript

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  library(ggpubr)
  library(argparse)
  library(Biostrings)
})

# 创建命令行参数解析器
create_parser <- function() {
  parser <- ArgumentParser(description = '线粒体突变分析工具 - 生成所有分析图表')
  
  parser$add_argument('--coverage', type = 'character', 
                     help = '覆盖度文件路径 (barcode_coverage.tsv.gz)')
  parser$add_argument('--variants', type = 'character',
                     help = '变异稀疏矩阵文件路径 (variant_sparse_matrix.tsv)')
  parser$add_argument('--somatic', type = 'character',
                     help = 'somatic突变文件路径 (snv.somatic.tsv)')
  parser$add_argument('--celltype', type = 'character',
                     help = '细胞类型文件路径 (human-mix-info.csv)')
  parser$add_argument('--bam', type = 'character',
                     help = 'BAM文件路径 (可选)')
  parser$add_argument('--barcodes', type = 'character',
                     help = 'barcode列表文件路径 (可选)')
  parser$add_argument('--reference', type = 'character',
                     help = '线粒体参考基因组文件路径')
  parser$add_argument('--output-dir', type = 'character', default = '.',
                     help = '输出目录路径')
  parser$add_argument('--prefix', type = 'character', default = 'analysis',
                     help = '输出文件前缀')
  parser$add_argument('--sample-cells', type = 'integer', default = 0,
                     help = '只处理指定数量的细胞进行测试 (0表示处理所有)')
  
  return(parser)
}

# A. Coverage分析函数
analyze_coverage <- function(coverage_file, celltype_file, output_dir, prefix) {
  message("A. 正在分析覆盖度...")
  
  # 读取覆盖度数据
  df <- read_tsv(coverage_file, col_names = c("barcode", "position", "coverage"))
  df <- df %>%
    mutate(
      position = as.integer(position),
      coverage = as.numeric(coverage)
    )
  
  message("读取了 ", nrow(df), " 行数据")
  message("包含 ", length(unique(df$barcode)), " 个细胞")
  
  # 读取细胞类型信息
  if (!is.null(celltype_file) && file.exists(celltype_file)) {
    celltype_df <- read_csv(celltype_file)
    
    if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
      celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
      message("加载了 ", length(celltype_mapping), " 个细胞的类型信息")
      
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
  
  # 计算每个细胞类型的平均覆盖度
  avg_coverage <- df %>%
    group_by(celltype, position) %>%
    summarise(avg_coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
    complete(celltype, position = 1:16569, fill = list(avg_coverage = 0))
  
  # 计算统计信息
  celltype_stats <- df %>%
    group_by(barcode, celltype) %>%
    summarise(total_coverage = sum(coverage), .groups = "drop") %>%
    group_by(celltype) %>%
    summarise(
      n_cells = n(),
      mean_coverage = mean(total_coverage / 16569),
      .groups = "drop"
    )
  
  # 绘制覆盖度图
  mt_regions <- data.frame(
    name = c("D-loop", "12S rRNA", "16S rRNA", "ND1", "ND2", "CO1", "CO2", 
             "ATP8", "ATP6", "CO3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"),
    start = c(1, 648, 1671, 3307, 4470, 5904, 7586, 8366, 8527, 9207, 
              10059, 10470, 10760, 12337, 14149, 14747),
    end = c(576, 1601, 3229, 4262, 5511, 7445, 8269, 8572, 9207, 9990, 
            10404, 10766, 12137, 14148, 14673, 15887),
    type = c("D-loop", "rRNA", "rRNA", "Protein", "Protein", "Protein", 
             "Protein", "Protein", "Protein", "Protein", "Protein", 
             "Protein", "Protein", "Protein", "Protein", "Protein")
  )
  
  plot_data <- avg_coverage %>%
    left_join(celltype_stats, by = "celltype") %>%
    mutate(celltype_label = paste0(celltype, " (n=", n_cells, ")"))
  
  max_coverage <- max(plot_data$avg_coverage, na.rm = TRUE)
  
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
    labs(title = "A. Mitochondrial Coverage by Cell Type",
         x = "Position (bp)", 
         y = "Coverage",
         color = "Cell Type",
         fill = "Genomic Region") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
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
  output_file <- file.path(output_dir, paste0(prefix, "_A_coverage.pdf"))
  ggsave(output_file, p, width = 14, height = 8)
  message("覆盖度图已保存至: ", output_file)
  
  return(list(plot = p, stats = celltype_stats))
}

# B. mtDNA百分比分析函数
analyze_mtdna_percent <- function(bam_file, barcodes_file, celltype_file, output_dir, prefix) {
  message("B. 正在分析mtDNA百分比...")
  
  # 使用samtools提取信息
  extract_bam_stats <- function(bam_file, barcode_file) {
    target_barcodes <- fread(barcode_file, header = FALSE)$V1
    
    cmd <- sprintf("samtools view %s | awk '{barcode=\"\"; for(i=12;i<=NF;i++) if($i~/^CB:Z:/) {barcode=substr($i,6); break} if(barcode!=\"\") print barcode\"\\t\"$3}'", bam_file)
    
    data <- fread(cmd = cmd, sep = "\t", header = FALSE, 
                  col.names = c("barcode", "chr"))
    
    filtered <- data[barcode %in% target_barcodes, ]
    
    total <- filtered[, .(total_reads = .N), by = barcode]
    mt <- filtered[chr %in% c("chrM", "MT", "M", "chrMT"), 
                   .(mt_reads = .N), by = barcode]
    
    result <- merge(total, mt, by = "barcode", all.x = TRUE)
    result[is.na(mt_reads), mt_reads := 0]
    result[, mt_ratio := mt_reads / total_reads * 100]
    
    return(result[order(-total_reads)])
  }
  
  stats <- extract_bam_stats(bam_file, barcodes_file)
  
  # 添加细胞类型信息
  if (!is.null(celltype_file) && file.exists(celltype_file)) {
    celltype_df <- read_csv(celltype_file)
    
    if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
      celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
      
      stats <- stats %>% filter(barcode %in% names(celltype_mapping))
      stats$celltype <- celltype_mapping[stats$barcode]
    } else {
      stats$celltype <- "Unknown"
    }
  } else {
    stats$celltype <- "Unknown"
  }
  
  # 绘制小提琴图
  p <- ggplot(stats, aes(x = celltype, y = mt_ratio, fill = celltype)) +
    geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
                 color = "black", linewidth = 0.4, fatten = 1.5) +
    geom_point(position = position_jitter(width = 0.15, height = 0), 
               alpha = 0.4, size = 1.2, color = "black", shape = 16) +
    labs(
      title = "B. mtDNA Content by Cell Type",
      x = "Cell Type",
      y = "mtDNA Content (%)"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text = element_text(family = "Arial", color = "black", size = 10),
      axis.title = element_text(family = "Arial", color = "black", size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  output_file <- file.path(output_dir, paste0(prefix, "_B_mtdna_percent.pdf"))
  ggsave(output_file, p, width = 11, height = 5.5)
  message("mtDNA百分比图已保存至: ", output_file)
  
  return(list(plot = p, stats = stats))
}

# C. 突变负担分析函数
analyze_mutation_burden <- function(variant_file, celltype_file, output_dir, prefix) {
  message("C. 正在分析突变负担...")
  
  # 读取变异数据
  variants <- fread(variant_file)
  
  # 计算每个细胞中VAF不为0的mutation数目
  cell_mutations <- variants %>%
    filter(vaf > 0) %>%
    group_by(cell) %>%
    summarise(
      mutation_count = n(),
      mean_vaf = mean(vaf, na.rm = TRUE),
      max_vaf = max(vaf, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 添加细胞类型信息
  if (!is.null(celltype_file) && file.exists(celltype_file)) {
    celltype_df <- fread(celltype_file)
    
    if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
      celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
      
      cell_mutations$celltype <- celltype_mapping[cell_mutations$cell]
      cell_mutations$celltype[is.na(cell_mutations$celltype)] <- "Unknown"
    } else {
      cell_mutations$celltype <- "Unknown"
    }
  } else {
    cell_mutations$celltype <- "Unknown"
  }
  
  # 绘制小提琴图
  p <- ggplot(cell_mutations, aes(x = celltype, y = mutation_count, fill = celltype)) +
    geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
                 color = "black", linewidth = 0.4, fatten = 1.5) +
    geom_point(position = position_jitter(width = 0.15, height = 0), 
               alpha = 0.4, size = 1.2, color = "black", shape = 16) +
    labs(
      title = "C. Mutation Burden by Cell Type",
      x = "Cell Type",
      y = "Number of Mutations per Cell (VAF > 0)"
    ) +
    theme_classic() +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  output_file <- file.path(output_dir, paste0(prefix, "_C_mutation_burden.pdf"))
  ggsave(output_file, p, width = 11, height = 5.5)
  message("突变负担图已保存至: ", output_file)
  
  return(list(plot = p, stats = cell_mutations))
}

# D. 标准化突变负担分析函数
analyze_normalized_burden <- function(variant_file, mtdna_stats, celltype_file, output_dir, prefix) {
  message("D. 正在分析标准化突变负担...")
  
  # 读取变异数据
  variants <- fread(variant_file)
  
  # 计算每个细胞中VAF不为0的mutation数目
  cell_mutations <- variants %>%
    filter(vaf > 0) %>%
    group_by(cell) %>%
    summarise(
      mutation_count = n(),
      .groups = "drop"
    )
  
  # 合并数据
  merged_data <- cell_mutations %>%
    rename(barcode = cell) %>%
    inner_join(mtdna_stats, by = "barcode") %>%
    distinct(barcode, .keep_all = TRUE)
  
  # 计算标准化突变负担
  merged_data <- merged_data %>%
    mutate(
      mutation_burden_norm = ifelse(mt_ratio > 0, mutation_count / mt_ratio, NA)
    )
  
  # 添加细胞类型信息
  if (!is.null(celltype_file) && file.exists(celltype_file)) {
    celltype_df <- fread(celltype_file)
    
    if ("barcode" %in% colnames(celltype_df) && "Annotation" %in% colnames(celltype_df)) {
      celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
      merged_data$celltype <- celltype_mapping[merged_data$barcode]
      merged_data$celltype[is.na(merged_data$celltype)] <- "Unknown"
    }
  }
  
  # 绘制小提琴图
  p <- ggplot(merged_data, aes(x = celltype, y = mutation_burden_norm, fill = celltype)) +
    geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
                 color = "black", linewidth = 0.4, fatten = 1.5) +
    geom_point(position = position_jitter(width = 0.15, height = 0), 
               alpha = 0.4, size = 1.2, color = "black", shape = 16) +
    labs(
      title = "D. Normalized Mutation Burden by Cell Type",
      x = "Cell Type",
      y = "Normalized Mutation Burden"
    ) +
    theme_classic() +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  output_file <- file.path(output_dir, paste0(prefix, "_D_normalized_burden.pdf"))
  ggsave(output_file, p, width = 11, height = 5.5)
  message("标准化突变负担图已保存至: ", output_file)
  
  return(list(plot = p, data = merged_data))
}

# E. 突变特征谱分析函数
analyze_mutation_spectrum <- function(somatic_file, reference_file, output_dir, prefix) {
  message("E. 正在分析突变特征谱...")
  
  # 读取somatic突变
  somatic_snv <- fread(somatic_file)
  colnames(somatic_snv) <- c("position","ref","alt","ref_fw","ref_rev","alt_fw",
                             "alt_rev","strand_score","mean_vaf","var_vaf","lis",
                             "pct_conf","pct_vaf_pos")
  
  mutation_list <- unique(paste0(somatic_snv$position, somatic_snv$ref, ">", somatic_snv$alt))
  
  # 简单的反向互补函数
  reverse_complement <- function(s) {
    chartr("ATGC", "TACG", s)
  }
  
  # 读取参考基因组
  fasta_seqs <- readDNAStringSet(reference_file)
  seq_names <- names(fasta_seqs)
  chrM_index <- grep("chrM|MT|mitochondria", seq_names, ignore.case = TRUE)
  chrM_seq <- as.character(fasta_seqs[[chrM_index]])
  
  ref_all <- data.frame(
    pos = 1:nchar(chrM_seq),
    ref = strsplit(chrM_seq, "")[[1]]
  )
  ref_all$ref <- toupper(ref_all$ref)
  
  # 构建三核苷酸上下文
  l <- as.character(ref_all$ref)
  ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))
  ref_all <- ref_all[!grepl("N", ref_all$three), ]
  
  # 创建所有可能的突变
  ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt, ]
  
  # 添加元数据
  ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, 
                                 ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))
  
  # 确定链
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")
  
  # 转换为C/T作为参考碱基
  ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
  ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", 
                                    ref_all_long$three, ref_all_long$rc3)
  ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", 
                                      ref_all_long$change, ref_all_long$change_rc)
  
  # 注释突变
  ref_all_long$called <- ref_all_long$variant %in% mutation_list
  
  # 计算期望/观察比例
  total <- dim(ref_all_long)[1]
  total_called <- sum(ref_all_long$called)
  
  prop_df <- ref_all_long %>% 
    group_by(three_plot, group_change, strand) %>%
    summarize(observed_prop_called = sum(called)/total_called, 
              expected_prop = n()/total, n = n()) %>%
    mutate(fc_called = observed_prop_called/expected_prop)
  
  prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
  
  # 绘制突变特征图
  p <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, size = 4)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_hline(yintercept = 1, linetype = 2, color = "black") +
    labs(title = "E. Mutation Spectrum",
         x = "Change in Nucleotide Context", 
         y = "Substitution Rate (Observed / Expected)") +
    facet_wrap(~ group_change, scales = "free_x", nrow = 1)
  
  output_file <- file.path(output_dir, paste0(prefix, "_E_mutation_spectrum.pdf"))
  ggsave(output_file, p, width = 8, height = 3)
  message("突变特征谱图已保存至: ", output_file)
  
  return(list(plot = p, data = prop_df))
}

# F. 突变频率谱分析函数
analyze_mutation_frequency <- function(variant_file, somatic_file, celltype_file, output_dir, prefix) {
  message("F. 正在分析突变频率谱...")
  
  # 读取数据
  variants <- fread(variant_file)
  somatic_snv <- fread(somatic_file)
  
  # 处理somatic突变列表
  colnames(somatic_snv) <- c("position","ref","alt","ref_fw","ref_rev","alt_fw",
                             "alt_rev","strand_score","mean_vaf","var_vaf","lis",
                             "pct_conf","pct_vaf_pos")
  
  mutation_list <- unique(paste0(somatic_snv$position, somatic_snv$ref, 
                                 ">", somatic_snv$alt))
  
  # 过滤变异
  variants <- variants[which(variants$vaf > 0), ]
  variants$snv <- paste0(variants$pos, paste(variants$ref_base, 
                                             variants$alt_base, sep = ">"))
  variants <- variants[which(variants$snv %in% mutation_list), ]
  
  # 添加细胞类型信息
  if (!is.null(celltype_file) && file.exists(celltype_file)) {
    celltype_df <- fread(celltype_file)
    celltype_mapping <- setNames(celltype_df$Annotation, celltype_df$barcode)
    
    variants_with_celltype <- variants %>%
      mutate(
        vaf_percent = vaf * 100,
        celltype = celltype_mapping[cell]
      ) %>%
      filter(!is.na(celltype) & celltype != "")
  } else {
    variants_with_celltype <- variants %>%
      mutate(
        vaf_percent = vaf * 100,
        celltype = "All"
      )
  }
  
  # 创建频率分组
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
  
  # 计算频率分布
  freq_dist <- variants_with_celltype %>%
    filter(!is.na(vaf_group)) %>%
    group_by(celltype, vaf_group) %>%
    summarise(
      mutation_count = n(),
      .groups = "drop"
    ) %>%
    complete(celltype, vaf_group, fill = list(mutation_count = 0)) %>%
    mutate(vaf_group = factor(vaf_group, levels = labels))
  
  # 绘制分面图
  p <- ggplot(freq_dist, aes(x = vaf_group, y = mutation_count, fill = vaf_group)) +
    geom_bar(stat = "identity", alpha = 0.8, color = "black") +
    facet_wrap(~ celltype, scales = "free_y", ncol = 3) +
    labs(
      title = "F. Mutation Frequency Distribution by Cell Type",
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
  
  output_file <- file.path(output_dir, paste0(prefix, "_F_mutation_frequency.pdf"))
  ggsave(output_file, p, width = 10, height = 10)
  message("突变频率谱图已保存至: ", output_file)
  
  return(list(plot = p, data = freq_dist))
}

# G. VAF与置信度分析函数
analyze_vaf_correlation <- function(somatic_file, output_dir, prefix) {
  message("G. 正在分析VAF与置信度相关性...")
  
  # 读取somatic突变
  somatic_snv <- fread(somatic_file)
  colnames(somatic_snv) <- c("position","ref","alt","ref_fw","ref_rev","alt_fw",
                             "alt_rev","strand_score","mean_vaf","var_vaf","lis",
                             "pct_conf","pct_vaf_pos")
  
  # 绘制VAF vs 置信度图
  p1 <- ggplot(somatic_snv, aes(x = mean_vaf, y = pct_conf)) +
    geom_point(alpha = 0.6, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top",
             size = 5, color = "black", fontface = "bold") +
    labs(
      title = "G. Mean VAF vs Confidence Percentage",
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
  
  output_file <- file.path(output_dir, paste0(prefix, "_G_vaf_correlation.pdf"))
  ggsave(output_file, p1, width = 8, height = 6)
  message("VAF相关性图已保存至: ", output_file)
  
  return(list(plot = p1, data = somatic_snv))
}

# H. LIS分布分析函数
analyze_lis_distribution <- function(somatic_file, output_dir, prefix) {
  message("H. 正在分析LIS分布...")
  
  # 读取somatic突变
  somatic_snv <- fread(somatic_file)
  colnames(somatic_snv) <- c("position","ref","alt","ref_fw","ref_rev","alt_fw",
                             "alt_rev","strand_score","mean_vaf","var_vaf","lis",
                             "pct_conf","pct_vaf_pos")
  
  # 绘制LIS密度分布图
  p <- ggplot(somatic_snv, aes(x = lis)) +
    geom_density(fill = "steelblue", alpha = 0.6, color = "steelblue", linewidth = 1) +
    geom_vline(aes(xintercept = median(lis, na.rm = TRUE)), 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = mean(lis, na.rm = TRUE)), 
               color = "darkgreen", linetype = "dashed", linewidth = 1) +
    annotate("text",
             x = quantile(somatic_snv$lis, 0.98, na.rm = TRUE),
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
      title = "H. Distribution of LIS Values",
      x = "LIS (Lineage Informative Score)",
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
  
  output_file <- file.path(output_dir, paste0(prefix, "_H_lis_distribution.pdf"))
  ggsave(output_file, p, width = 10, height = 7)
  message("LIS分布图已保存至: ", output_file)
  
  return(list(plot = p, stats = summary(somatic_snv$lis)))
}

# 创建组合图函数
create_composite_plot <- function(all_plots, output_dir, prefix) {
  message("正在创建组合图...")
  
  # 提取所有图形
  plots <- list(
    all_plots$coverage$plot,
    all_plots$mtdna_percent$plot,
    all_plots$mutation_burden$plot,
    all_plots$normalized_burden$plot,
    all_plots$mutation_spectrum$plot,
    all_plots$mutation_frequency$plot,
    all_plots$vaf_correlation$plot,
    all_plots$lis_distribution$plot
  )
  
  # 创建布局
  composite_plot <- wrap_plots(plots, ncol = 2, nrow = 4) +
    plot_annotation(
      title = "Comprehensive Mitochondrial Mutation Analysis",
      subtitle = "A-H: Complete analysis pipeline",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5)
      )
    )
  
  # 保存组合图
  output_file <- file.path(output_dir, paste0(prefix, "_composite_plot.pdf"))
  ggsave(output_file, composite_plot, width = 20, height = 24)
  message("组合图已保存至: ", output_file)
  
  return(composite_plot)
}

# 主函数
main <- function() {
  # 解析命令行参数
  parser <- create_parser()
  args <- parser$parse_args()
  
  # 创建输出目录
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  # 存储所有分析结果
  all_results <- list()
  
  # 执行所有分析（检查文件是否存在）
  tryCatch({
    # A. Coverage分析
    if (!is.null(args$coverage) && file.exists(args$coverage)) {
      all_results$coverage <- analyze_coverage(args$coverage, args$celltype, 
                                               args$output_dir, args$prefix)
    } else {
      message("跳过Coverage分析：文件不存在")
    }
    
    # B. mtDNA百分比分析
    if (!is.null(args$bam) && !is.null(args$barcodes) && 
        file.exists(args$bam) && file.exists(args$barcodes)) {
      all_results$mtdna_percent <- analyze_mtdna_percent(args$bam, args$barcodes, 
                                                         args$celltype, args$output_dir, 
                                                         args$prefix)
    } else {
      message("跳过mtDNA百分比分析：BAM或barcodes文件不存在")
    }
    
    # C. 突变负担分析
    if (!is.null(args$variants) && file.exists(args$variants)) {
      all_results$mutation_burden <- analyze_mutation_burden(args$variants, 
                                                            args$celltype, 
                                                            args$output_dir, 
                                                            args$prefix)
    } else {
      message("跳过突变负担分析：变异文件不存在")
    }
    
    # D. 标准化突变负担分析（需要mtDNA统计数据）
    if (!is.null(args$variants) && file.exists(args$variants) && 
        !is.null(all_results$mtdna_percent)) {
      all_results$normalized_burden <- analyze_normalized_burden(args$variants, 
                                                                all_results$mtdna_percent$stats,
                                                                args$celltype, 
                                                                args$output_dir, 
                                                                args$prefix)
    } else {
      message("跳过标准化突变负担分析：所需数据不完整")
    }
    
    # E. 突变特征谱分析
    if (!is.null(args$somatic) && !is.null(args$reference) && 
        file.exists(args$somatic) && file.exists(args$reference)) {
      all_results$mutation_spectrum <- analyze_mutation_spectrum(args$somatic, 
                                                                args$reference, 
                                                                args$output_dir, 
                                                                args$prefix)
    } else {
      message("跳过突变特征谱分析：somatic或参考基因组文件不存在")
    }
    
    # F. 突变频率谱分析
    if (!is.null(args$variants) && !is.null(args$somatic) && 
        file.exists(args$variants) && file.exists(args$somatic)) {
      all_results$mutation_frequency <- analyze_mutation_frequency(args$variants, 
                                                                  args$somatic, 
                                                                  args$celltype, 
                                                                  args$output_dir, 
                                                                  args$prefix)
    } else {
      message("跳过突变频率谱分析：所需文件不存在")
    }
    
    # G. VAF与置信度分析
    if (!is.null(args$somatic) && file.exists(args$somatic)) {
      all_results$vaf_correlation <- analyze_vaf_correlation(args$somatic, 
                                                            args$output_dir, 
                                                            args$prefix)
    } else {
      message("跳过VAF相关性分析：somatic文件不存在")
    }
    
    # H. LIS分布分析
    if (!is.null(args$somatic) && file.exists(args$somatic)) {
      all_results$lis_distribution <- analyze_lis_distribution(args$somatic, 
                                                              args$output_dir, 
                                                              args$prefix)
    } else {
      message("跳过LIS分布分析：somatic文件不存在")
    }
    
    # 创建组合图（如果至少有2个分析成功）
    if (length(all_results) >= 2) {
      all_results$composite_plot <- create_composite_plot(all_results, 
                                                         args$output_dir, 
                                                         args$prefix)
    }
    
    # 保存分析结果汇总
    saveRDS(all_results, file.path(args$output_dir, 
                                   paste0(args$prefix, "_analysis_results.rds")))
    
    message("\n=== 分析完成 ===")
    message("输出目录: ", args$output_dir)
    message("文件前缀: ", args$prefix)
    message("完成的分析: ", paste(names(all_results), collapse = ", "))
    
  }, error = function(e) {
    message("分析过程中出现错误: ", e$message)
    return(1)
  })
  
  return(0)
}

# 运行主函数
if (interactive()) {
  # 在RStudio中交互式运行
  message("在交互式模式下运行，请提供参数")
} else {
  # 作为脚本运行
  exit_code <- main()
  quit(save = "no", status = exit_code)
}