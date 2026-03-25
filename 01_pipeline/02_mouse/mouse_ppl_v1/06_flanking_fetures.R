# 继续前面的代码...

# 1. 首先将突变分为6种类型（基于ref和alt）
ref_all_long <- ref_all_long %>%
  mutate(
    mutation_class = group_change
  )


ref_all_long <- ref_all_long %>% filter(!is.na(mutation_class))
ref_all_long <- ref_all_long [which(ref_all_long$called=="TRUE"),]

# 提取每个突变位置上下5bp的序列
extract_context <- function(position, flank = 5) {
  start_pos <- max(1, position - flank)
  end_pos <- min(nchar(chrM_seq), position + flank)
  context_seq <- substr(chrM_seq, start_pos, end_pos)
  return(context_seq)
}

# 提取中心位置（突变位置）上下5bp的序列
ref_all_long$context_seq <- sapply(ref_all_long$pos, extract_context, flank = 5)

# 2. 定义分析特征的函数
analyze_sequence_features <- function(seq) {
  # 将序列转为大写
  seq <- toupper(seq)
  # 序列长度
  seq_len <- nchar(seq)
  # 中心位置索引（第6个字符，因为是5bp flanking）
  center_pos <- 6
  
  # 提取不同区域
  upstream <- substr(seq, 1, 5)  # 上游5bp
  downstream <- substr(seq, 7, 11)  # 下游5bp
  left_flank <- substr(seq, 1, center_pos - 1)  # 左侧flanking
  right_flank <- substr(seq, center_pos + 1, seq_len)  # 右侧flanking
  
  # 计算GC含量
  calc_gc_content <- function(s) {
    if (nchar(s) == 0) return(0)
    chars <- strsplit(s, "")[[1]]
    gc_count <- sum(chars %in% c("G", "C"))
    return(gc_count / nchar(s))
  }
  
  # 计算AT含量
  calc_at_content <- function(s) {
    if (nchar(s) == 0) return(0)
    chars <- strsplit(s, "")[[1]]
    at_count <- sum(chars %in% c("A", "T"))
    return(at_count / nchar(s))
  }
  
  # 检测重复序列模式（简单重复检测）
  detect_repeats <- function(s, min_repeat = 2) {
    if (nchar(s) < 2) return(FALSE)
    
    # 检查单核苷酸重复
    chars <- strsplit(s, "")[[1]]
    max_single_repeat <- max(rle(chars)$lengths)
    
    # 检查二核苷酸重复
    di_nuc_repeat <- FALSE
    if (nchar(s) >= 4) {
      di_nucs <- sapply(1:(nchar(s)-1), function(i) substr(s, i, i+1))
      max_di_repeat <- max(rle(di_nucs)$lengths)
      di_nuc_repeat <- max_di_repeat >= min_repeat
    }
    
    return(max_single_repeat >= 3 | di_nuc_repeat)  # 单核苷酸重复≥3次或二核苷酸重复≥2次
  }
  
  # 计算嘌呤含量
  calc_purine_content <- function(s) {
    if (nchar(s) == 0) return(0)
    chars <- strsplit(s, "")[[1]]
    purine_count <- sum(chars %in% c("A", "G"))
    return(purine_count / nchar(s))
  }
  
  # 计算嘧啶含量
  calc_pyrimidine_content <- function(s) {
    if (nchar(s) == 0) return(0)
    chars <- strsplit(s, "")[[1]]
    pyrimidine_count <- sum(chars %in% c("C", "T"))
    return(pyrimidine_count / nchar(s))
  }
  
  # 计算序列熵（复杂度）
  calc_entropy <- function(s) {
    if (nchar(s) == 0) return(0)
    chars <- strsplit(s, "")[[1]]
    freq <- table(chars) / length(chars)
    entropy <- -sum(freq * log2(freq))
    return(entropy)
  }
  
  # 返回所有特征
  return(data.frame(
    # 长度信息
    total_length = seq_len,
    
    # 含量特征
    gc_content_total = calc_gc_content(seq),
    gc_content_upstream = calc_gc_content(upstream),
    gc_content_downstream = calc_gc_content(downstream),

    at_content_total = calc_at_content(seq),
    at_content_upstream = calc_at_content(upstream),
    at_content_downstream = calc_at_content(downstream),
    
    purine_content_total = calc_purine_content(seq),
    purine_content_upstream = calc_purine_content(upstream),
    purine_content_downstream = calc_purine_content(downstream),
    
    pyrimidine_content_total = calc_pyrimidine_content(seq),
    pyrimidine_content_upstream = calc_pyrimidine_content(upstream),
    pyrimidine_content_downstream = calc_pyrimidine_content(downstream),
    
    # 重复特征
    has_repeat_upstream = detect_repeats(upstream),
    has_repeat_downstream = detect_repeats(downstream),
    has_repeat_total = detect_repeats(seq),
    
    # 复杂度特征
    entropy_total = calc_entropy(seq),
    entropy_upstream = calc_entropy(upstream),
    entropy_downstream = calc_entropy(downstream),
    
    # 序列本身
    upstream_seq = upstream,
    downstream_seq = downstream,
    center_base = substr(seq, center_pos, center_pos),
    
    stringsAsFactors = FALSE
  ))
}

# 3. 为每个突变计算特征
features_list <- lapply(ref_all_long$context_seq, analyze_sequence_features)
features_df <- do.call(rbind, features_list)

# 合并到原数据框
ref_all_long_features <- cbind(ref_all_long, features_df)

# 4. 按突变类型分组分析
mutation_type_analysis <- ref_all_long_features %>%
  group_by(mutation_class) %>%
  summarize(
    # 基本统计
    count = n(),
    
    # GC含量统计
    mean_gc_total = mean(gc_content_total, na.rm = TRUE),
    mean_gc_upstream = mean(gc_content_upstream, na.rm = TRUE),
    mean_gc_downstream = mean(gc_content_downstream, na.rm = TRUE),
    sd_gc_total = sd(gc_content_total, na.rm = TRUE),
    
    # AT含量统计
    mean_at_total = mean(at_content_total, na.rm = TRUE),
    mean_at_upstream = mean(at_content_upstream, na.rm = TRUE),
    mean_at_downstream = mean(at_content_downstream, na.rm = TRUE),
    
    # 嘌呤/嘧啶含量
    mean_purine_total = mean(purine_content_total, na.rm = TRUE),
    mean_pyrimidine_total = mean(pyrimidine_content_total, na.rm = TRUE),
    
    # 重复序列比例
    prop_repeat_upstream = mean(has_repeat_upstream, na.rm = TRUE),
    prop_repeat_downstream = mean(has_repeat_downstream, na.rm = TRUE),
    prop_repeat_total = mean(has_repeat_total, na.rm = TRUE),
    
    # 序列复杂度
    mean_entropy_total = mean(entropy_total, na.rm = TRUE),
    mean_entropy_upstream = mean(entropy_upstream, na.rm = TRUE),
    mean_entropy_downstream = mean(entropy_downstream, na.rm = TRUE),
    
    # 是否存在突变
    called_count = sum(called),
    called_prop = sum(called) / n()
  ) %>%
  arrange(desc(count))

# 显示结果
print(mutation_type_analysis)




# 5. Visualization
library(ggplot2)
library(tidyr)
library(patchwork)  # 用于多图组合
library(dplyr)

# 创建一个函数来生成并保存所有特征的图表
data=mutation_type_analysis_clean
raw_data=ref_all_long_features_clean
output_dir = "./"


  # 1. Mutation count distribution (放在最前面)
  count_data <- data %>%
    select(mutation_class, count) %>%
    filter(!is.na(mutation_class))
  
  p1 <- ggplot(count_data, aes(x = reorder(mutation_class, -count), y = count, fill = mutation_class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::comma(count)), vjust = -0.5, size = 4) +
    labs(title = "Mutation Count by Type",
         x = "Mutation Type",
         y = "Count",
         fill = "Mutation Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1)))
  
  # 2. GC content violin plot (每个mutation的分布)
  # 首先需要准备原始数据
  if("gc_content_total" %in% colnames(raw_data)) {
    gc_raw_data <- raw_data %>%
      filter(!is.na(mutation_class)) %>%
      select(mutation_class, gc_content_total)
    
    p2 <- ggplot(gc_raw_data, aes(x = mutation_class, y = gc_content_total, fill = mutation_class)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      labs(title = "GC Content Distribution by Mutation Type",
           x = "Mutation Type",
           y = "GC Content",
           fill = "Mutation Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none")
  }
  
  # 3. AT content violin plot
  if("at_content_total" %in% colnames(raw_data)) {
    at_raw_data <- raw_data %>%
      filter(!is.na(mutation_class)) %>%
      select(mutation_class, at_content_total)
    
    p3 <- ggplot(at_raw_data, aes(x = mutation_class, y = at_content_total, fill = mutation_class)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      labs(title = "AT Content Distribution by Mutation Type",
           x = "Mutation Type",
           y = "AT Content",
           fill = "Mutation Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none")
  }
  
  # 4. Sequence entropy violin plot
  if("entropy_total" %in% colnames(raw_data)) {
    entropy_raw_data <- raw_data %>%
      filter(!is.na(mutation_class)) %>%
      select(mutation_class, entropy_total)
    
    p4 <- ggplot(entropy_raw_data, aes(x = mutation_class, y = entropy_total, fill = mutation_class)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      labs(title = "Sequence Entropy Distribution by Mutation Type",
           x = "Mutation Type",
           y = "Entropy (bits)",
           fill = "Mutation Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none")
  }
  
  # 5. Repeat sequence analysis - 识别频繁出现的repeat序列
  # 假设raw_data中有上游和下游序列信息

    # 提取所有repeat序列（这里以重复长度≥3的序列为例）

find_repeats <- function(seq) {
  if(is.na(seq) || nchar(seq) < 3) return(NA)
  repeats <- c()
  # 1. 单核苷酸重复 (A{3,}, C{3,}, G{3,}, T{3,})
  single_repeats <- gregexpr("([ACGT])\\1{2,}", seq, perl = TRUE)
  matches <- regmatches(seq, single_repeats)[[1]]
  if(length(matches) > 0) {
    repeats <- c(repeats, matches)
  }
    # 2. 二核苷酸重复 (AT){2,}, (CG){2,} 等
  if(nchar(seq) >= 4) {
    di_repeats <- gregexpr("([ACGT]{2})\\1{1,}", seq, perl = TRUE)
    di_matches <- regmatches(seq, di_repeats)[[1]]
    if(length(di_matches) > 0) {
      repeats <- c(repeats, di_matches)
    }
  }
  
  # 3. 三核苷酸重复 (CAG){2,}, (GAA){2,} 等
  if(nchar(seq) >= 6) {
    tri_repeats <- gregexpr("([ACGT]{3})\\1{1,}", seq, perl = TRUE)
    tri_matches <- regmatches(seq, tri_repeats)[[1]]
    if(length(tri_matches) > 0) {
      repeats <- c(repeats, tri_matches)
    }
  }
  if(length(repeats) == 0) return(NA)
  return(paste(repeats, collapse = ";"))
}


    # 分析上游和下游的repeat序列
    raw_data <- raw_data %>%
      filter(!is.na(mutation_class)) %>%
      mutate(
        upstream_repeats = sapply(upstream_seq, find_repeats),
        downstream_repeats = sapply(downstream_seq, find_repeats)
      )
    
    # 合并所有repeat序列
    all_repeats <- raw_data %>%
      select(mutation_class, upstream_repeats, downstream_repeats) %>%
      pivot_longer(cols = c(upstream_repeats, downstream_repeats), 
                   names_to = "region", values_to = "repeat_seq") %>%
      filter(!is.na(repeat_seq)) %>%
      separate_rows(repeat_seq, sep = ";") %>%
      filter(nchar(repeat_seq) >= 3)  # 只保留长度≥3的重复
    
    # 统计每个mutation类型中最常见的repeat序列
    top_repeats <- all_repeats %>%
      group_by(mutation_class, repeat_seq) %>%
      summarize(
        count = n(),
        .groups = 'drop'
      ) %>%
      group_by(mutation_class) %>%
      slice_max(order_by = count, n = 10000) %>%  # 每个类型取前10
      ungroup()
    
    # 保存repeat序列统计
    write.csv(top_repeats, paste0(output_dir, "top_repeat_sequences.csv"), row.names = FALSE)
    
    # 创建热图数据
    # 选择每个类型中top 5的repeat序列用于热图
    heatmap_repeats <- top_repeats %>%
      group_by(mutation_class) %>%
      slice_max(order_by = count, n = 1000) %>%
      ungroup() %>%
      arrange(mutation_class, desc(count))
    
    # 创建热图
    p5 <- ggplot(heatmap_repeats, aes(x = mutation_class, y = reorder(repeat_seq, count), fill = count)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = count), color = "white", size = 3) +
      scale_fill_gradientn(
        colors = c("#2E86AB", "#A23B72", "#F18F01"),
        name = "Count",
        trans = "log10"  # 使用对数转换以更好显示差异
      ) +
      labs(title = "Top Repeat Sequences by Mutation Type",
           x = "Mutation Type",
           y = "Repeat Sequence") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid = element_blank()
      )
    
    # 6. Repeat sequence length distribution
    repeat_length_data <- all_repeats %>%
      mutate(repeat_length = nchar(repeat_seq))
    
    p6 <- ggplot(repeat_length_data, aes(x = mutation_class, y = repeat_length, fill = mutation_class)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      labs(title = "Repeat Sequence Length Distribution",
           x = "Mutation Type",
           y = "Repeat Length (bp)",
           fill = "Mutation Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none")
    
    # 7. Repeat type analysis (单核苷酸重复 vs 多核苷酸重复)
    repeat_type_data <- all_repeats %>%
      mutate(
        repeat_type = ifelse(
          nchar(repeat_seq) == nchar(gsub("(.)\\1+", "\\1", repeat_seq)) * length(strsplit(repeat_seq, "")[[1]]) / nchar(repeat_seq),
          "Homopolymer",
          "Other"
        )
      ) %>%
      group_by(mutation_class, repeat_type) %>%
      summarize(count = n(), .groups = 'drop') %>%
      group_by(mutation_class) %>%
      mutate(prop = count / sum(count))
    
    p7 <- ggplot(repeat_type_data, aes(x = mutation_class, y = prop, fill = repeat_type)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(aes(label = paste0(round(prop * 100, 1), "%")),
                position = position_stack(vjust = 0.5), size = 3) +
      labs(title = "Repeat Sequence Type Distribution",
           x = "Mutation Type",
           y = "Proportion",
           fill = "Repeat Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "bottom")
  
  # 保存单个PDF
  pdf(paste0(output_dir, "mutation_feature_analysis.pdf"), width = 5, height = 4)
  print(p1)
  if(exists("p2")) print(p2)
  if(exists("p3")) print(p3)
  if(exists("p4")) print(p4)
  if(exists("p5")) print(p5)
  if(exists("p6")) print(p6)
  if(exists("p7")) print(p7) 
  dev.off()
   
   
# 6. 补充分析：重复序列的具体模式
# 识别最常见的重复序列模式
cat("\n=== Repeat Sequence Pattern Analysis ===\n")

# 假设我们已经有了all_repeats数据
if(exists("all_repeats")) {
  # 分析最常见的重复序列
  top_global_repeats <- all_repeats %>%
    group_by(repeat_seq) %>%
    summarize(
      total_count = n(),
      mutation_types = paste(unique(mutation_class), collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(desc(total_count)) %>%
    slice_head(n = 20)
  
  print("Top 20 most frequent repeat sequences:")
  print(top_global_repeats)
  
  # 保存结果
  write.csv(top_global_repeats, "top_global_repeat_sequences.csv", row.names = FALSE)
  
  # 绘制全局重复序列分布
  p_top_repeats <- ggplot(top_global_repeats %>% slice_head(n = 15), 
                         aes(x = reorder(repeat_seq, total_count), y = total_count)) +
    geom_bar(stat = "identity", fill = "#2E86AB") +
    geom_text(aes(label = total_count), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(title = "Top 15 Most Frequent Repeat Sequences",
         x = "Repeat Sequence",
         y = "Count") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  ggsave("top_global_repeats.png", p_top_repeats, width = 10, height = 8, dpi = 300)
}

# 7. 保存所有结果
write.csv(mutation_type_analysis_clean, "mutation_type_sequence_features_clean.csv", row.names = FALSE)

# 8. 生成摘要报告
cat("\n=== Mutation Feature Analysis Summary ===\n")
cat("\n1. Mutation Statistics:\n")
for(mtype in mutation_type_analysis_clean$mutation_class) {
  count <- mutation_type_analysis_clean %>% 
    filter(mutation_class == mtype) %>% 
    pull(count)
  cat(sprintf("%s: %s mutations\n", mtype, scales::comma(count)))
}

cat("\n2. Analysis Complete. Files Generated:\n")
cat("- mutation_feature_analysis.pdf (comprehensive PDF report)\n")
cat("- Individual PNG files for each plot\n")
cat("- mutation_type_sequence_features_clean.csv\n")
if(exists("all_repeats")) {
  cat("- top_repeat_sequences.csv\n")
  cat("- top_global_repeat_sequences.csv\n")
  cat("- top_global_repeats.png\n")
}

cat("\nAnalysis completed successfully!\n")