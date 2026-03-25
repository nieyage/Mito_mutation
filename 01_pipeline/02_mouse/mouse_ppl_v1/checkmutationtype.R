mut_data_high_con <- read.csv("mut_data_high_con_rescue_CA.csv")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(stringr)
library(Biostrings)
library(data.table)
library(tidyr)

# 创建突变类型的基础函数
reverse_complement <- function(s) {
  chartr("ATGC", "TACG", s)
}

# 加载参考基因组
unshifted_chrM_ref <- "/md01/nieyg/ref/mito_ref/mm10/mm10.chrM.fasta"
fasta_seqs <- readDNAStringSet(unshifted_chrM_ref)

# 查找chrM序列
seq_names <- names(fasta_seqs)
chrM_index <- grep("chrM|MT|mitochondria", seq_names, ignore.case = TRUE)
chrM_seq <- as.character(fasta_seqs[[chrM_index]])

# 构建参考基因组数据
ref_all <- data.frame(
  pos = 1:nchar(chrM_seq),
  ref = strsplit(chrM_seq, "")[[1]]
)
ref_all$ref <- toupper(ref_all$ref)

# 创建三碱基上下文
l <- as.character(ref_all$ref)
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))
ref_all <- ref_all[!grepl("N", ref_all$three),]

# 创建所有可能的突变
ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# 添加元数据
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# 确定链方向
table(ref_all$ref)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C", "T"), "L", "H")

# 转换为C/T作为参考等位基因
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# 计算总突变位点数和不同类型突变的可能发生概率
total_sites <- nrow(ref_all_long)
total_by_type <- ref_all_long %>% 
  group_by(three_plot, group_change, strand) %>%
  summarize(n_possible = n(), .groups = 'drop') %>%
  mutate(possible_prob = n_possible / total_sites)

plot_mutation_signature <- function(cells, title, plot_type = "standard") {
  # 提取这些细胞的突变
  subset_data <- mut_data_high_con %>% 
    filter(barcode %in% cells) %>%
    mutate(mutation_type = paste0(Ref, ">", VarAllele),
           variant = paste0(Position, mutation_type))
  
  # 获取突变列表
  called_variants <- unique(subset_data$variant)
  
  # 标记哪些位点发生了突变
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  # 计算总观测突变数
  total_observed <- length(called_variants)
  
  if (total_observed == 0) {
    message(paste("Cells", paste(cells, collapse=","), "have no mutations detected"))
    return(NULL)
  }
  
  # 根据plot_type选择不同的计算方法
  if (plot_type == "standard") {
    # 第一部分：使用原始方法计算fc_called
    total_possible <- nrow(ref_all_long)
    
    prop_df <- ref_all_long %>% 
      group_by(three_plot, group_change, strand) %>%
      summarize(
        observed_prop_called = sum(called) / total_observed,
        expected_prop = n() / total_possible,
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        fc_called = observed_prop_called / expected_prop,
        change_plot = paste0(group_change, "_", three_plot)
      )
    
    y_var <- "fc_called"
    y_label <- "Observed proportion / Expected proportion"
    plot_subtitle <- paste0("Total mutations: ", total_observed)
    
  } else if (plot_type == "observed_vs_genomic") {
    # 第二部分：观测概率 vs 基因组可能发生概率
    total_possible <- nrow(ref_all_long)
    
    # 首先，统计观察到的突变中各种类型的数量
    subset_data$three_plot<- ref_all_long[match(subset_data$variant,ref_all_long$variant),]$three_plot

    observed_counts <- subset_data %>%
      mutate(
        change = paste0(Ref, ">", VarAllele),
        strand = ifelse(Ref %in% c("C", "T"), "L", "H"),
        group_change = ifelse(strand == "L", 
                             paste0(Ref, VarAllele),
                             reverse_complement(paste0(Ref, VarAllele)))
      ) %>%
      group_by(three_plot, group_change, strand) %>%
      summarize(observed_count = n(), .groups = 'drop')
    
    # 总观察突变数（不同位点）
    total_observed_variants <- nrow(subset_data)
    
    # 合并基因组信息
    genomic_counts <- ref_all_long %>%
      group_by(three_plot, group_change, strand) %>%
      summarize(
        genomic_count = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        genomic_prob = genomic_count / total_possible
      )
    prop_df<- as.data.frame(observed_counts)
    prop_df$observed_prob<- prop_df$observed_count / total_observed_variants
    prop_df$label<- paste0(prop_df$three_plot,prop_df$group_change,prop_df$strand)
    genomic_counts<- as.data.frame(genomic_counts)
    genomic_counts$label<- paste0(genomic_counts$three_plot,genomic_counts$group_change,genomic_counts$strand)
    prop_df2 <-merge(prop_df,genomic_counts,by="label")
    prop_df2$ratio <- prop_df2$observed_prob/prop_df2$genomic_prob
    prop_df2$change_plot <- paste0(prop_df2$group_change.x, "_", prop_df2$three_plot.x)
    prop_df2$group_change <- prop_df2$group_change.x
    prop_df2$strand <- prop_df2$strand.x

    # 将观察计数与基因组计数合并
    prop_df <- prop_df2
    y_var <- "ratio"
    y_label <- "Observed probability / Genomic probability"
    plot_subtitle <- paste0("Observed mutations: ", total_observed_variants, 
                           " | Unique mutation types: ", nrow(observed_counts))
    
  } else if (plot_type == "absolute_counts") {
    # 第三部分：绝对计数
    subset_data$three_plot<- ref_all_long[match(subset_data$variant,ref_all_long$variant),]$three_plot
    absolute_counts <- subset_data %>%
      mutate(
        change = paste0(Ref, ">", VarAllele),
        strand = ifelse(Ref %in% c("C", "T"), "L", "H"),
        group_change = ifelse(strand == "L", 
                             paste0(Ref, VarAllele),
                             reverse_complement(paste0(Ref, VarAllele)))
      ) %>%
      group_by(three_plot, group_change, strand) %>%
      summarize(absolute_count = n(), .groups = 'drop')

    absolute_counts$change_plot<- paste0(absolute_counts$group_change, "_", absolute_counts$three_plot)
    
    # 合并数据
    prop_df <- absolute_counts
    y_var <- "absolute_count"
    y_label <- "Mutation count"
    plot_subtitle <- paste0("Total mutations: ", total_observed, 
                           " | Unique variants: ", length(called_variants))
  }
  
  # 绘图
  if (plot_type == "absolute_counts") {
    # 对于绝对计数，使用连续颜色
    p <- ggplot(prop_df, aes(x = change_plot, y = !!sym(y_var), fill = !!sym(y_var))) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                         name = "Count") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80")
      ) +
      labs(
        x = "Mutation type (trinucleotide context)",
        y = y_label,
        title = title,
        subtitle = plot_subtitle
      ) +
      facet_wrap(~ group_change, scales = "free_x", nrow = 1)
    
  } else {
    # 对于比例图，使用strand作为颜色
    p <- ggplot(prop_df, aes(x = change_plot, y = !!sym(y_var), fill = strand)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
      scale_fill_manual(values = c("firebrick", "dodgerblue3"), 
                        name = "Strand",
                        labels = c("H" = "Heavy", "L" = "Light")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80")
      )
    
    if (plot_type == "standard") {
      p <- p + 
        geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(
          x = "Mutation type (trinucleotide context)",
          y = y_label,
          title = title,
          subtitle = plot_subtitle
        )
    } else {
      p <- p + 
        geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(
          x = "Mutation type (trinucleotide context)",
          y = y_label,
          title = title,
          subtitle = plot_subtitle
        )
    }
    
    p <- p + facet_wrap(~ group_change, scales = "free_x", nrow = 1) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  }
  
  return(list(plot = p, data = prop_df, plot_type = plot_type))
}


set.seed(1234)

all_cells <- unique(mut_data_high_con$barcode)
selected_cells_1 <- sample(all_cells, min(1, length(all_cells)))
selected_cells_50 <- sample(all_cells, min(50, length(all_cells)))
selected_cells_100 <- sample(all_cells, min(100, length(all_cells)))


library(gridExtra)
library(grid)
pdf("./50cell_for_signature.pdf", width = 14, height = 12)

  result_50_standard <- plot_mutation_signature(selected_cells_50, 
                                               "50 Cells - Standard Method", 
                                               plot_type = "standard")
  result_50_obs_gen <- plot_mutation_signature(selected_cells_50, 
                                              "50 Cells - Observed vs Genomic", 
                                              plot_type = "observed_vs_genomic")
  result_50_abs <- plot_mutation_signature(selected_cells_50, 
                                          "50 Cells - Absolute Counts", 
                                          plot_type = "absolute_counts")
  all_plots <- list(
    result_50_standard$plot, result_50_obs_gen$plot, result_50_abs$plot  )
# 创建标题
  grid.arrange(
    arrangeGrob(
      textGrob("Mutation Signature Analysis Comparison", 
               gp = gpar(fontsize = 18, fontface = "bold")),
      arrangeGrob(
        grobs = all_plots,
        ncol = 1,
        nrow = 3
      ),
      ncol = 1,
      heights = c(0.05, 0.95)
    )
  )
  dev.off()

library(gridExtra)
library(grid)
pdf("./allcells_for_signature.pdf", width = 14, height = 12)

  result_50_standard <- plot_mutation_signature(all_cells, 
                                               "All Cells - Standard Method", 
                                               plot_type = "standard")
  result_50_obs_gen <- plot_mutation_signature(all_cells, 
                                              "All Cells - Observed vs Genomic", 
                                              plot_type = "observed_vs_genomic")
  result_50_abs <- plot_mutation_signature(all_cells, 
                                          "All Cells - Absolute Counts", 
                                          plot_type = "absolute_counts")
  all_plots <- list(
    result_50_standard$plot, result_50_obs_gen$plot, result_50_abs$plot  )
# 创建标题
  grid.arrange(
    arrangeGrob(
      textGrob("Mutation Signature Analysis Comparison", 
               gp = gpar(fontsize = 18, fontface = "bold")),
      arrangeGrob(
        grobs = all_plots,
        ncol = 1,
        nrow = 3
      ),
      ncol = 1,
      heights = c(0.05, 0.95)
    )
  )
  dev.off()

library(gridExtra)

cell_counts <- c(1, 20, 50,100,8900)
all_cells <- unique(mut_data_high_con$barcode)

for (n in cell_counts) {
  if (n > length(all_cells)) next
  
  pdf(paste0("mutation_", n, "cells.pdf"), 14, 36)
  
  plots <- list()
  
  for (r in 1:3) {
    set.seed(123 + r * n)
    cells <- sample(all_cells, n)
    
    p_std <- plot_mutation_signature(cells, "", "standard")
    p_obs <- plot_mutation_signature(cells, "", "observed_vs_genomic")  
    p_abs <- plot_mutation_signature(cells, "", "absolute_counts")
    
    if (!is.null(p_std)) plots <- c(plots, list(p_std$plot))
    if (!is.null(p_obs)) plots <- c(plots, list(p_obs$plot))
    if (!is.null(p_abs)) plots <- c(plots, list(p_abs$plot))
  }
  
  if (length(plots) > 0) {
    grid.arrange(grobs = plots, ncol = 1)
  }
  
  dev.off()
}
