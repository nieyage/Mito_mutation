plot_mutation_signature_4in1 <- function(mutation_list, title_prefix = "Mutation Signature Analysis") {
  # mutation_list: 数据框，必须包含以下列：Position, Ref, VarAllele
  
  # 验证输入数据
  required_cols <- c("Position", "Ref", "VarAllele")
  if (!all(required_cols %in% colnames(mutation_list))) {
    stop(paste("Input data must contain columns:", paste(required_cols, collapse = ", ")))
  }
  
  # 去重复，每个突变位点只计数一次（不同细胞中相同突变视为一个突变）
  unique_mutations <- mutation_list %>%
    distinct(Position, Ref, VarAllele) %>%
    mutate(
      mutation_type = paste0(Ref, ">", VarAllele),
      variant = paste0(Position, mutation_type)
    )
  
  # 获取突变列表
  called_variants <- unique(unique_mutations$variant)
  
  if (length(called_variants) == 0) {
    message("No mutations detected in the input list")
    return(NULL)
  }
  
  # 标记哪些位点发生了突变
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  # 计算总观测突变数（去重后）
  total_observed <- length(called_variants)
  
  message(paste("Total unique mutations:", total_observed))
  
  # 总可能突变数
  total_possible <- nrow(ref_all_long)
  
  # =====================================================================
  # 数据准备：图1和图2的数据（基于去重的唯一突变）
  # =====================================================================
  
  # 计算每个突变类型的观测比例和预期比例（唯一突变）
  prop_df_unique <- ref_all_long %>% 
    group_by(three_plot, group_change, strand) %>%
    summarize(
      observed_count_unique = sum(called),  # 去重计数
      observed_prop_unique = sum(called) / total_observed,
      expected_prop = n() / total_possible,
      n_possible = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      ratio_unique = observed_prop_unique / expected_prop,
      change_plot = paste0(group_change, "_", three_plot)
    )
  
  # =====================================================================
  # 数据准备：图3和图4的数据（包含所有突变）
  # =====================================================================
  
  # 这次考虑所有突变（包括不同细胞中的相同突变）
  all_mutations <- mutation_list %>%
    mutate(
      mutation_type = paste0(Ref, ">", VarAllele),
      variant = paste0(Position, mutation_type)
    )
  
  # 获取突变类型的三碱基上下文
  all_mutations$three_plot <- ref_all_long[match(all_mutations$variant, ref_all_long$variant),]$three_plot
  all_mutations$group_change <- ref_all_long[match(all_mutations$variant, ref_all_long$variant),]$group_change
  all_mutations$strand <- ref_all_long[match(all_mutations$variant, ref_all_long$variant),]$strand
  
  # 统计所有突变的类型（包括重复）
  observed_counts_all <- all_mutations %>%
    group_by(three_plot, group_change, strand) %>%
    summarize(
      observed_count_all = n(),  # 包括重复
      .groups = 'drop'
    )
  
  total_observed_all <- nrow(all_mutations)
  message(paste("Total mutations (including duplicates):", total_observed_all))
  
  # 计算观测概率（包含重复）
  observed_counts_all$observed_prob_all <- observed_counts_all$observed_count_all / total_observed_all
  
  # 获取基因组概率
  genomic_prob <- ref_all_long %>%
    group_by(three_plot, group_change, strand) %>%
    summarize(
      genomic_count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      genomic_prob = genomic_count / total_possible
    )
  
  # 合并数据（包含重复）
  prop_df_all <- observed_counts_all %>%
    left_join(genomic_prob, by = c("three_plot", "group_change", "strand")) %>%
    mutate(
      ratio_all = observed_prob_all / genomic_prob,
      change_plot = paste0(group_change, "_", three_plot)
    ) %>%
    filter(!is.na(ratio_all))  # 移除没有匹配的行
  
  # =====================================================================
  # 创建基础绘图函数
  # =====================================================================
  
  create_mutation_plot <- function(data, y_var, y_label, plot_title, show_hline = FALSE) {
    # 获取突变类型分组
    mutation_groups <- unique(data$group_change)
    n_groups <- length(mutation_groups)
    
    p <- ggplot(data, aes(x = change_plot, y = !!sym(y_var), fill = strand)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
      scale_fill_manual(
        values = c("firebrick", "dodgerblue3"), 
        name = "Strand",
        labels = c("H" = "Heavy", "L" = "Light")
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(
          angle = 90, 
          hjust = 1, 
          vjust = 0.5, 
          size = 5.5,
          margin = margin(t = 2)
        ),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 5)),
        axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 5)),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 5)),
        legend.position = "none",  # 移除单个图例
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80", linewidth = 0.3),
        plot.margin = margin(5, 5, 5, 5),
        strip.text = element_text(size = 7, face = "bold", margin = margin(1, 0, 1, 0)),
        strip.background = element_rect(fill = "grey95", color = "grey80", linewidth = 0.3)
      ) +
      labs(
        x = "Mutation type (trinucleotide context)",
        y = y_label,
        title = plot_title
      ) +
      facet_wrap(
        ~ group_change, 
        scales = "free_x", 
        nrow = 1,
        labeller = labeller(group_change = function(x) gsub("(\\w{2})", "\\1>", x))
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),
        limits = if (show_hline) c(0, max(data[[y_var]], na.rm = TRUE) * 1.1) else NULL
      )
    
    if (show_hline) {
      p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.3)
    }
    
    return(p)
  }
  
  # =====================================================================
  # 创建四个子图
  # =====================================================================
  
  # 图1: 观测比例/预期比例（唯一突变）
  p1 <- create_mutation_plot(
    data = prop_df_unique,
    y_var = "ratio_unique",
    y_label = "Observed/Expected\n(unique)",
    plot_title = "A. Ratio of Proportions (Unique Mutations)",
    show_hline = TRUE
  )
  
  # 图2: 观测绝对计数（唯一突变）
  p2 <- create_mutation_plot(
    data = prop_df_unique,
    y_var = "observed_count_unique",
    y_label = "Mutation count\n(unique)",
    plot_title = "B. Absolute Counts (Unique Mutations)",
    show_hline = FALSE
  )
  
  # 图3: 观测概率/基因组概率（所有突变）
  p3 <- create_mutation_plot(
    data = prop_df_all,
    y_var = "ratio_all",
    y_label = "Observed/Genomic\n(all)",
    plot_title = "C. Probability Ratio (All Mutations)",
    show_hline = TRUE
  )
  
  # 图4: 观测绝对计数（所有突变）
  p4 <- create_mutation_plot(
    data = prop_df_all,
    y_var = "observed_count_all",
    y_label = "Mutation count\n(all)",
    plot_title = "D. Absolute Counts (All Mutations)",
    show_hline = FALSE
  )
  
  # =====================================================================
  # 创建共享图例
  # =====================================================================
  
  # 创建一个简单的图例图
  legend_data <- data.frame(
    strand = factor(c("H", "L"), levels = c("H", "L")),
    value = c(1, 1)
  )
  
  legend_plot <- ggplot(legend_data, aes(x = strand, y = value, fill = strand)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("firebrick", "dodgerblue3"), 
      name = "Strand",
      labels = c("H" = "Heavy (A/G rich)", "L" = "Light (C/T rich)")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold", margin = margin(b = 5)),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      legend.key = element_rect(fill = "white", color = "grey50"),
      legend.box.margin = margin(5, 0, 5, 0)
    ) +
    guides(fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    ))
  
  # 提取图例
  shared_legend <- cowplot::get_legend(legend_plot)
  
  # =====================================================================
  # 组合四个子图
  # =====================================================================
  
  # 使用 patchwork 包组合图形
  if (!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
  # 组合四个子图（四行）
  combined_plot <- p1 / p2 / p3 / p4 / shared_legend +
    plot_layout(
      heights = c(1, 1, 1, 1, 0.1),  # 四行图 + 图例
      guides = 'collect'
    ) &
    theme(
      legend.position = 'none'  # 移除各个子图的图例
    )
  
  # 添加整体标题和注释
  combined_plot <- combined_plot +
    plot_annotation(
      title = paste(title_prefix, "- Comprehensive Mutation Signature Analysis"),
      subtitle = paste0(
        "Unique mutations: ", total_observed, 
        " | Total occurrences: ", total_observed_all,
        " | Average occurrences per site: ", round(total_observed_all / total_observed, 2)
      ),
      caption = paste0(
        "Plots A & B: Based on unique mutations (positions counted once)\n",
        "Plots C & D: Based on all mutation occurrences (including duplicates)\n",
        "Expected/Genomic probability: Proportion of possible mutation sites in reference genome"
      ),
      theme = theme(
        plot.title = element_text(
          size = 16, 
          face = "bold", 
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.subtitle = element_text(
          size = 11, 
          hjust = 0.5,
          margin = margin(b = 10)
        ),
        plot.caption = element_text(
          size = 9, 
          hjust = 0,
          color = "gray40",
          margin = margin(t = 10)
        ),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  # =====================================================================
  # 汇总统计信息
  # =====================================================================
  
  summary_stats <- list(
    unique_mutations = total_observed,
    total_mutations = total_observed_all,
    mutation_types = length(unique(all_mutations$mutation_type)),
    avg_occurrences_per_site = round(total_observed_all / total_observed, 2),
    data_unique = prop_df_unique,
    data_all = prop_df_all
  )
  
  # 打印汇总信息
  cat("=== Mutation Signature Analysis Summary ===\n")
  cat(paste("• Unique mutations (distinct positions):", total_observed, "\n"))
  cat(paste("• Total mutation occurrences:", total_observed_all, "\n"))
  cat(paste("• Mutation types:", length(unique(all_mutations$mutation_type)), "\n"))
  cat(paste("• Average occurrences per site:", round(total_observed_all / total_observed, 2), "\n"))
  cat(paste("• Most frequent mutation type:\n"))
  
  # 找出最常见的突变类型
  if (nrow(prop_df_all) > 0) {
    most_frequent <- prop_df_all[which.max(prop_df_all$observed_count_all), ]
    cat(paste("   - Type:", most_frequent$group_change, "\n"))
    cat(paste("   - Context:", most_frequent$three_plot, "\n"))
    cat(paste("   - Count:", most_frequent$observed_count_all, "\n"))
    cat(paste("   - Strand:", ifelse(most_frequent$strand == "L", "Light", "Heavy"), "\n"))
  }
  
  return(list(
    combined_plot = combined_plot,
    individual_plots = list(
      plot_ratio_unique = p1,
      plot_count_unique = p2,
      plot_ratio_all = p3,
      plot_count_all = p4
    ),
    stats = summary_stats,
    data = list(
      unique_mutations = unique_mutations,
      all_mutations = all_mutations,
      summary_unique = prop_df_unique,
      summary_all = prop_df_all
    )
  ))
}

# 使用示例函数
plot_and_save_mutation_signature <- function(mutation_list, output_prefix = "mutation_signature", 
                                             width = 16, height = 20, ...) {
  # 运行分析函数
  result <- plot_mutation_signature_4in1(mutation_list, ...)
  
  if (is.null(result)) {
    message("No mutations to plot")
    return(NULL)
  }
  
  # 保存组合图
  combined_filename <- paste0(output_prefix, "_4in1.pdf")
  ggsave(
    combined_filename,
    result$combined_plot,
    width = width,
    height = height,
    device = "pdf"
  )
  message(paste("Saved combined plot to:", combined_filename))
  
  # 保存单个图（可选）
  single_dir <- paste0(output_prefix, "_individual_plots")
  if (!dir.exists(single_dir)) dir.create(single_dir)
  
  ggsave(
    file.path(single_dir, "1_ratio_unique.pdf"),
    result$individual_plots$plot_ratio_unique,
    width = width,
    height = height/4 + 2,
    device = "pdf"
  )
  
  ggsave(
    file.path(single_dir, "2_count_unique.pdf"),
    result$individual_plots$plot_count_unique,
    width = width,
    height = height/4 + 2,
    device = "pdf"
  )
  
  ggsave(
    file.path(single_dir, "3_ratio_all.pdf"),
    result$individual_plots$plot_ratio_all,
    width = width,
    height = height/4 + 2,
    device = "pdf"
  )
  
  ggsave(
    file.path(single_dir, "4_count_all.pdf"),
    result$individual_plots$plot_count_all,
    width = width,
    height = height/4 + 2,
    device = "pdf"
  )
  
  # 保存统计数据
  stats_filename <- paste0(output_prefix, "_statistics.txt")
  sink(stats_filename)
  cat("=== Mutation Signature Analysis Statistics ===\n\n")
  cat(paste("Analysis Date:", Sys.Date(), "\n"))
  cat(paste("Output Prefix:", output_prefix, "\n\n"))
  
  cat("Summary Statistics:\n")
  cat(paste("Unique mutations (distinct positions):", result$stats$unique_mutations, "\n"))
  cat(paste("Total mutation occurrences:", result$stats$total_mutations, "\n"))
  cat(paste("Mutation types:", result$stats$mutation_types, "\n"))
  cat(paste("Average occurrences per site:", result$stats$avg_occurrences_per_site, "\n\n"))
  
  cat("Top 10 Most Frequent Mutation Types (all occurrences):\n")
  if (nrow(result$data$summary_all) > 0) {
    top10 <- result$data$summary_all %>%
      arrange(desc(observed_count_all)) %>%
      head(10) %>%
      select(group_change, three_plot, strand, observed_count_all, ratio_all)
    
    print(top10)
  }
  sink()
  message(paste("Saved statistics to:", stats_filename))
  
  # 返回结果
  return(result)
}

# 使用示例：
# mutation_data <- data.frame(
#   Position = c(100, 200, 300, 100, 200, 400, 500, 100, 200),
#   Ref = c("A", "C", "G", "A", "C", "T", "A", "A", "C"),
#   VarAllele = c("G", "T", "A", "G", "T", "C", "G", "G", "T")
# )
# 
# result <- plot_and_save_mutation_signature(
#   mutation_data,
#   output_prefix = "my_mutation_analysis",
#   title_prefix = "My Sample"
# )
# 
# # 显示组合图
# print(result$combined_plot)