python /md01/nieyg/pipeline/mito_mutation/01_human/SNV_summary_v3.py -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/output/Donor1 -o ./Donor1_snv_v3.tsv -v 

python3 /md01/nieyg/pipeline/mito_mutation/01_human/filter_mutations.py \
  -i ./Donor1_snv_v3.tsv \
  -o ./ \
  --germline-output germline.tsv \
  --somatic-output somatic.tsv \
  --strand-min 0.3 \
  --strand-max 0.7 \
  --min-depth 10 \
  --alt-ratio-threshold 0.90 \
  --no-compress


# 3. variant_sparse_matrix.tsv

python3 /md01/nieyg/pipeline/mito_mutation/01_human/create_variant_matrix_v2.py \
-s somatic.tsv \
-d /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_testppl_output/output/Donor1 \
-o variant_sparse_matrix.tsv.gz --workers 20




library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

somatic_snv<- "./somatic.tsv"
somatic_snv <- fread(somatic_snv)


head(somatic_snv)

# 安装并加载ggpubr包
# install.packages("ggpubr")
library(ggpubr)

p1 <- ggplot(somatic_snv, aes(x = Mean_vaf, y = Pct_conf)) +
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

p2 <- ggplot(somatic_snv, aes(x = Mean_vaf, y = Pct_vaf_pos)) +
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


# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggpubr)

# 定义6种突变类型（CA、CG、CT、TA、TC、TG），考虑互补链
somatic_snv <- somatic_snv %>%
  mutate(
    Mutation_Type = case_when(
      (Ref == "C" & Alt == "A") | (Ref == "G" & Alt == "T") ~ "CA",
      (Ref == "C" & Alt == "G") | (Ref == "G" & Alt == "C") ~ "CG",
      (Ref == "C" & Alt == "T") | (Ref == "G" & Alt == "A") ~ "CT",
      (Ref == "T" & Alt == "A") | (Ref == "A" & Alt == "T") ~ "TA",
      (Ref == "T" & Alt == "C") | (Ref == "A" & Alt == "G") ~ "TC",
      (Ref == "T" & Alt == "G") | (Ref == "A" & Alt == "C") ~ "TG",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Mutation_Type != "Other")  # 过滤掉其他类型

# 设置颜色
colors <- c(CA = "#E41A1C", CG = "#377EB8", CT = "#4DAF4A", 
            TA = "#984EA3", TC = "#FF7F00", TG = "#FFFF33")

# 创建因子
somatic_snv$Mutation_Type <- factor(somatic_snv$Mutation_Type, 
                                   levels = names(colors))

# 绘制散点图（整体趋势 + 颜色分组）
plot_scatter <- function(data, x, y, title) {
  ggplot(data, aes(x = !!sym(x), y = !!sym(y), color = Mutation_Type)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", 
                linetype = "dashed", size = 0.8, aes(group = 1)) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = colors) +
    stat_cor(label.x.npc = "left", label.y.npc = "top") +
    labs(title = title, x = paste0(x, " (log10)"), y = paste0(y, " (log10)")) +
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# 生成两个主要图形
p1 <- plot_scatter(somatic_snv, "Mean_vaf", "Pct_conf", 
                   "Mean VAF vs Confidence Percentage")
p2 <- plot_scatter(somatic_snv, "Mean_vaf", "Pct_vaf_pos", 
                   "Mean VAF vs VAF-positive Cells Percentage")

# 保存图形
ggsave("plot1_mean_vaf_vs_pct_conf.pdf", p1, width = 9, height = 7)
ggsave("plot2_mean_vaf_vs_pct_vaf_pos.pdf", p2, width = 9, height = 7)

# 输出统计摘要
summary_stats <- somatic_snv %>%
  group_by(Mutation_Type) %>%
  summarise(n = n(), 
            mean_VAF = mean(Mean_vaf),
            mean_Conf = mean(Pct_conf),
            mean_VAF_pos = mean(Pct_vaf_pos),
            .groups = "drop")

print(summary_stats)
write.csv(summary_stats, "mutation_summary.csv", row.names = FALSE)


