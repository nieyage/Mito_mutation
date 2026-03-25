# check the human ppl 

# the overlap of 4 ppl mutations 

# 

# 2.sc_SNV filter
# define file path and an empty dataframe
library(dplyr)
arg_1 = "/md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/masked_SNVcalling_percell/"
files <- dir(arg_1, pattern = "snv$")
path <- arg_1
i <- 1
final <- data.frame("Chrom" = 0 , "Position" = 0, "Ref" = 0, "VarAllele" = 0)
final <- final[-length(final$Chrom),]
sc_germline <- final

final <- data.frame()
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
  if(nrow(z)>0){
  z$barcode<- files[i]
  k <- filter(z, Reads2/(Reads1 + Reads2) >= 0.90)
  n <- filter(z, Reads2/(Reads1 + Reads2) < 0.90) #select(z, c("Chrom", "Position", "Ref", "VarAllele"))
  k <- k #select(k, c("Chrom", "Position", "Ref", "VarAllele"))
  final <- merge(final, n, all = T)}
  # final <- rbind(final, n)
  #sc_germline <- rbind(sc_germline, k)
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

# 

csv_file="/md01/nieyg/project/mito_mutation/01_pipeline/04_germline_mutation/human-mix-info.csv"  # 请替换为实际的CSV文件路径
metadata<- read.csv(csv_file)
mut_data$celltype<- metadata[match(mut_data$barcode,metadata$barcode),]$Annotation
mut_data$donor<- metadata[match(mut_data$barcode,metadata$barcode),]$sample
mut_data_high_con<- mut_data[which(mut_data$mutation_id%in%high_con_mutation),]
mut_data_high_con[order(mut_data_high_con$barcode),]

high_con_mutation<- paste(SNV_filter$Position, SNV_filter$Ref, SNV_filter$VarAllele, sep = "_")

mutation_1<- high_con_mutation


# jinxu ppl 
/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.somatic.tsv
/md01/jinxu/Project/mgatk-speedup/44_GenoByCell/variant_sparse_matrix.tsv

tmp2<- read.table("/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.somatic.tsv",row.names=1)
mutation_2<- paste(rownames(tmp2), tmp2$V2, tmp2$V3, sep = "_")


# human ppl v5 min reads =1 

tmp3_1<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor1/Donor1_somatic.tsv",sep="\t")
tmp3_2<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor2/Donor2_somatic.tsv",sep="\t")
tmp3_3<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor3/Donor3_somatic.tsv",sep="\t")
tmp3_4<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor4/Donor4_somatic.tsv",sep="\t")


mutation_3<- unique(c(paste(tmp3_1$Position, tmp3_1$Ref, tmp3_1$Alt, sep = "_"),
	paste(tmp3_2$Position, tmp3_2$Ref, tmp3_2$Alt, sep = "_"),
	paste(tmp3_3$Position, tmp3_3$Ref, tmp3_3$Alt, sep = "_"),
	paste(tmp3_4$Position, tmp3_4$Ref, tmp3_4$Alt, sep = "_")
	))




# human ppl v5 min reads =3

tmp4_1<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor1/Donor1_somatic.tsv",sep="\t")
tmp4_2<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor2/Donor2_somatic.tsv",sep="\t")
tmp4_3<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor3/Donor3_somatic.tsv",sep="\t")
tmp4_4<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor4/Donor4_somatic.tsv",sep="\t")
mutation_4<- unique(c(paste(tmp4_1$Position, tmp4_1$Ref, tmp4_1$Alt, sep = "_"),
	paste(tmp4_2$Position, tmp4_2$Ref, tmp4_2$Alt, sep = "_"),
	paste(tmp4_3$Position, tmp4_3$Ref, tmp4_3$Alt, sep = "_"),
	paste(tmp4_4$Position, tmp4_4$Ref, tmp4_4$Alt, sep = "_")
	))




# 安装并加载所需的包
if (!require("UpSetR")) install.packages("UpSetR")
if (!require("ggplot2")) install.packages("ggplot2")

human_ppl_v1<- mutation_1
jinxu_ppl<- mutation_2
human_ppl_v5_minreads1<- mutation_3
human_ppl_v5_minreads3<- mutation_4
library(UpSetR)
library(ggplot2)
all_mutations <- unique(c(mutation_1, mutation_2, mutation_3, mutation_4))
binary_matrix <- data.frame(
  mutation = all_mutations,
  human_ppl_v1 = as.integer(all_mutations %in% human_ppl_v1),
  jinxu_ppl = as.integer(all_mutations %in% jinxu_ppl),
  human_ppl_v5_minreads1 = as.integer(all_mutations %in% human_ppl_v5_minreads1),
  human_ppl_v5_minreads3 = as.integer(all_mutations %in% human_ppl_v5_minreads3)
)


pdf("mutations_upset_plot_4method.pdf", width = 8, height = 4)
upset(
  binary_matrix,
  sets = c("human_ppl_v1", "jinxu_ppl", "human_ppl_v5_minreads1", "human_ppl_v5_minreads3"),
  sets.bar.color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  main.bar.color = "#1F78B4",
  matrix.color = "#333333",
  order.by = "freq",
  empty.intersections = "on",
  mb.ratio = c(0.7, 0.3),  # 主图与矩阵图的比例
  text.scale = c(1.5, 1.2, 1, 1, 1.5, 1.2)  # 调整文本大小
)
dev.off()


extract_intersections <- function(binary_matrix, set_names) {
  results <- list()
  
  # 生成所有可能的组合（1到n个集合）
  n <- length(set_names)
  all_combos <- list()
  
  # 1个集合的组合
  for(i in 1:n) {
    all_combos[[paste(set_names[i])]] <- c(i)
  }
  
  # 2个集合的组合
  if(n >= 2) {
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        combo_name <- paste(set_names[i], set_names[j], sep = "&")
        all_combos[[combo_name]] <- c(i, j)
      }
    }
  }
  
  # 3个集合的组合
  if(n >= 3) {
    for(i in 1:(n-2)) {
      for(j in (i+1):(n-1)) {
        for(k in (j+1):n) {
          combo_name <- paste(set_names[i], set_names[j], set_names[k], sep = "&")
          all_combos[[combo_name]] <- c(i, j, k)
        }
      }
    }
  }
  
  # 4个集合的组合
  if(n >= 4) {
    for(i in 1:(n-3)) {
      for(j in (i+1):(n-2)) {
        for(k in (j+1):(n-1)) {
          for(l in (k+1):n) {
            combo_name <- paste(set_names[i], set_names[j], set_names[k], set_names[l], sep = "&")
            all_combos[[combo_name]] <- c(i, j, k, l)
          }
        }
      }
    }
  }
  
  # 计算每个组合的交集
  for(combo_name in names(all_combos)) {
    indices <- all_combos[[combo_name]]
    
    # 获取对应的列名
    cols <- set_names[indices]
    
    # 筛选出在这些列中值都为1的行
    mask <- rowSums(binary_matrix[, cols, drop = FALSE]) == length(cols)
    mutations <- binary_matrix$mutation[mask]
    
    results[[combo_name]] <- list(
      set_names = cols,
      count = length(mutations),
      mutations = mutations
    )
  }
  
  return(results)
}

# 使用函数
set_names <- c("human_ppl_v1", "jinxu_ppl", "human_ppl_v5_minreads1", "human_ppl_v5_minreads3")
intersection_results <- extract_intersections(binary_matrix, set_names)

# 查看所有组合
cat("所有交集组合统计:\n")
for(name in names(intersection_results)) {
  cat(sprintf("%-50s: %2d个突变\n", name, intersection_results[[name]]$count))
}
所有交集组合统计:
human_ppl_v1                                      : 3165个突变
jinxu_ppl                                         : 1549个突变
human_ppl_v5_minreads1                            : 3044个突变
human_ppl_v5_minreads3                            : 2804个突变
human_ppl_v1&jinxu_ppl                            : 1240个突变
human_ppl_v1&human_ppl_v5_minreads1               : 2613个突变
human_ppl_v1&human_ppl_v5_minreads3               : 2467个突变
jinxu_ppl&human_ppl_v5_minreads1                  : 1200个突变
jinxu_ppl&human_ppl_v5_minreads3                  : 997个突变
human_ppl_v5_minreads1&human_ppl_v5_minreads3     : 2573个突变
human_ppl_v1&jinxu_ppl&human_ppl_v5_minreads1     : 1149个突变
human_ppl_v1&jinxu_ppl&human_ppl_v5_minreads3     : 973个突变
human_ppl_v1&human_ppl_v5_minreads1&human_ppl_v5_minreads3: 2264个突变
jinxu_ppl&human_ppl_v5_minreads1&human_ppl_v5_minreads3: 992个突变
human_ppl_v1&jinxu_ppl&human_ppl_v5_minreads1&human_ppl_v5_minreads3: 970个突变

intersection_results[["human_ppl_v1&jinxu_ppl&human_ppl_v5_minreads1&human_ppl_v5_minreads3"]]


overlap <- intersection_results[["human_ppl_v1&jinxu_ppl&human_ppl_v5_minreads1&human_ppl_v5_minreads3"]]$mutations

# overlap mutation in human ppl v5 


library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

tmp3_1<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor1/Donor1_somatic.tsv",sep="\t")
tmp3_2<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor2/Donor2_somatic.tsv",sep="\t")
tmp3_3<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor3/Donor3_somatic.tsv",sep="\t")
tmp3_4<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor4/Donor4_somatic.tsv",sep="\t")

somatic_snv_v5<- rbind(tmp3_1,tmp3_2,tmp3_3,tmp3_4)

somatic_snv_v5$name<- paste(somatic_snv_v5$Position, somatic_snv_v5$Ref, somatic_snv_v5$Alt, sep = "_")

somatic_snv_v5<- somatic_snv_v5[which(somatic_snv_v5$name%in% overlap),]
somatic_snv_v5<- somatic_snv_v5[order(somatic_snv_v5$name),]



head(somatic_snv_v5)

# 安装并加载ggpubr包
# install.packages("ggpubr")
library(ggpubr)

p1 <- ggplot(somatic_snv_v5, aes(x = Mean_vaf, y = Pct_conf)) +
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

p2 <- ggplot(somatic_snv_v5, aes(x = Mean_vaf, y = Pct_vaf_pos)) +
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
ggsave("pplv5_overlap_mean_vaf_vs_pct_conf.pdf", p1, width = 8, height = 6)
ggsave("pplv5_overlap_mean_vaf_vs_pct_vaf_pos.pdf", p2, width = 8, height = 6)





# overlap mutation in jinxu ppl
# G. max VAF vs.  population score  

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

somatic_snv<- "/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.somatic.tsv"
somatic_snv <- fread(somatic_snv)
colnames(somatic_snv)<- c("position","ref","alt","ref_fw","ref_rev","alt_fw","alt_rev","strand_score","Mean_vaf","var_vaf","lis","Pct_conf","Pct_vaf_pos")
somatic_snv$name<- paste(somatic_snv$position, somatic_snv$ref, somatic_snv$alt, sep = "_")

somatic_snv<- somatic_snv[which(somatic_snv$name%in% overlap),]

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
ggsave("jinxuppl_overlap_mean_vaf_vs_pct_conf.pdf", p1, width = 8, height = 6)
ggsave("jinxuppl_overlap_mean_vaf_vs_pct_vaf_pos.pdf", p2, width = 8, height = 6)


# compare 
somatic_snv_v5<- somatic_snv_v5[order(somatic_snv_v5$name),]
somatic_snv<- somatic_snv[order(somatic_snv$name),]

head(somatic_snv_v5)
head(somatic_snv)
> somatic_snv[1:2,]
   position    ref    alt ref_fw ref_rev alt_fw alt_rev strand_score Mean_vaf
      <int> <char> <char>  <int>   <int>  <int>   <int>        <num>    <num>
1:     1000      T      C 107062   96530     47      33       0.0875    4e-04
2:    10015      T      C 112859   85365     14       8       0.1364    1e-04
   var_vaf   lis Pct_conf Pct_vaf_pos      name
     <num> <num>    <num>       <num>    <char>
1:   8e-05 4e-04    2e-04      0.0046  1000_T_C
2:   1e-05 1e-04    1e-04      0.0015 10015_T_C


> somatic_snv_v5[1:4,]
     Position Ref Alt Ref_fw_total Ref_rev_total Alt_fw_total Alt_rev_total
56       1000   T   C          115           103            7             8
1164     1000   T   C           73            69           12            11
3077     1000   T   C           29            33            9             6
3594    10015   T   C           22            18            4             3
     Strand_score Mean_vaf  Var_vaf      Lis Pct_conf Pct_vaf_pos
56         0.4667 0.073738 0.001120 0.073655    2e-04      0.0018
1164       0.5217 0.143471 0.016381 0.141159    4e-04      0.0015
3077       0.6000 0.226275 0.054303 0.214621    3e-04      0.0013
3594       0.5714 0.162050 0.007735 0.160806    3e-04      0.0006
     N_cells_vaf_pos Conf_cells Total_cells total_reads alt_reads  alt_ratio
56                 8          1        4369         233        15 0.06437768
1164               7          2        4591         165        23 0.13939394
3077               4          1        3098          77        15 0.19480519
3594               2          1        3098          47         7 0.14893617
          name
56    1000_T_C
1164  1000_T_C
3077  1000_T_C
3594 10015_T_C





overlap_3 <- intersection_results[["human_ppl_v1&human_ppl_v5_minreads1&human_ppl_v5_minreads3"]]$mutations
overlap_3<- overlap_3[-which(overlap_3 %in% c(jinxu_ppl))]

9942_G_A
685_A_T
7270_T_C

Donor1_variant_sparse_matrix<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor1/Donor1_variant_sparse_matrix_withoutstrand.tsv.gz",sep="\t")
Donor2_variant_sparse_matrix<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor2/Donor2_variant_sparse_matrix_withoutstrand.tsv.gz",sep="\t")
Donor3_variant_sparse_matrix<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor3/Donor3_variant_sparse_matrix_withoutstrand.tsv.gz",sep="\t")
Donor4_variant_sparse_matrix<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor4/Donor4_variant_sparse_matrix_withoutstrand.tsv.gz",sep="\t")


685%in% Donor1_variant_sparse_matrix$pos
685%in% Donor2_variant_sparse_matrix$pos
685%in% Donor3_variant_sparse_matrix$pos
685%in% Donor4_variant_sparse_matrix$pos
Donor2_variant_sparse_matrix[Donor2_variant_sparse_matrix$pos==685&Donor2_variant_sparse_matrix$vaf>0,]


Donor1_variant_sparse_matrix[Donor1_variant_sparse_matrix$pos==9942&Donor1_variant_sparse_matrix$vaf>0,]
                      cell  pos ref_base alt_base ref_count alt_count    vaf
964023  AGTGTGGCATAATTGC-1 9942        G        A        22         2 0.0833
1412130 GGATTGCGTTACAAAC-1 9942        G        A        10         2 0.1667
1700592 TTCGCAACAGGCCAAA-1 9942        G        A        19         2 0.0952
1949418 CAAGGCTGTATCTGGA-1 9942        G        A        16         3 0.1579
1964832 CAGCTAAGTGCAATGC-1 9942        G        A        15         1 0.0625
1998963 CGCATATAGCTATTAG-1 9942        G        A        24         2 0.0769
2014377 CCGTTTGGTAATGGCC-1 9942        G        A        33         1 0.0294
3058125 CTGTTTAGTTAAGCGC-1 9942        G        A        21        11 0.3438
4526859 CAAACTGGTTAGCAGC-1 9942        G        A        18         1 0.0526

tmp3_1<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor1/Donor1_somatic.tsv",sep="\t")
tmp3_2<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor2/Donor2_somatic.tsv",sep="\t")
tmp3_3<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor3/Donor3_somatic.tsv",sep="\t")
tmp3_4<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results_minreads1/Donor4/Donor4_somatic.tsv",sep="\t")

somatic_snv_v5<- rbind(tmp3_1,tmp3_2,tmp3_3,tmp3_4)

somatic_snv<- "/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.tsv"
somatic_snv <- fread(somatic_snv)
colnames(somatic_snv)<- c("position","ref","alt","ref_fw","ref_rev","alt_fw","alt_rev","strand_score","Mean_vaf","var_vaf","lis","Pct_conf","Pct_vaf_pos")
somatic_snv$name<- paste(somatic_snv$position, somatic_snv$ref, somatic_snv$alt, sep = "_")


somatic_snv_v5[which(somatic_snv_v5$Position==9942),]

grep "9942" AGTGTGGCATAATTGC-1.counts
grep "9942" GGATTGCGTTACAAAC-1.counts
grep "9942" TTCGCAACAGGCCAAA-1.counts
grep "9942" CAAGGCTGTATCTGGA-1.counts
grep "9942" CAGCTAAGTGCAATGC-1.counts
grep "9942" CGCATATAGCTATTAG-1.counts

 0.1187  0.008174        0.117738        0.0002  0.0021  9     

somatic_snv[which(somatic_snv$position==9942),]



overlap_4 <- intersection_results[["human_ppl_v5_minreads1&human_ppl_v5_minreads3"]]$mutations
overlap_4<- overlap_4[-which(overlap_4 %in% c(jinxu_ppl,human_ppl_v1))]

example:
5448_C_A
5448%in% Donor1_variant_sparse_matrix$pos
5448%in% Donor2_variant_sparse_matrix$pos
5448%in% Donor3_variant_sparse_matrix$pos
5448%in% Donor4_variant_sparse_matrix$pos
Donor2_variant_sparse_matrix[Donor2_variant_sparse_matrix$pos==5448&Donor2_variant_sparse_matrix$vaf>0,]

head -n 5450 ./Donor2/TCTTCAAGTGCGCAAT-1.counts |tail -n 20

# not in jinxu ppl 
somatic_snv[which(somatic_snv$position==5448),]

TCTTCAAGTGCGCAAT-1
# CA 和 GT 突变被过滤

13352_T_C
7331%in% Donor1_variant_sparse_matrix$pos
7331%in% Donor2_variant_sparse_matrix$pos
7331%in% Donor3_variant_sparse_matrix$pos
7331%in% Donor4_variant_sparse_matrix$pos

Donor3_variant_sparse_matrix[Donor3_variant_sparse_matrix$pos==7331&Donor3_variant_sparse_matrix$vaf>0,]
head -n 7331 ./Donor3/GCGGTTGGTGACATAT-1.counts |tail -n 50

head -n 7331 /md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/masked_SNVcalling_percell/GCGGTTGGTGACATAT-1.counts |tail -n 10

  x <- paste0(path, "AAGCGTTTCGTTTCCA-1.snv")
  y <- read.table(file = x, header = T, colClasses = c("character"))
y[y$Position==16482,]

                Reads2/(Reads1 + Reads2) >= 0.10



Donor2_variant_sparse_matrix[Donor2_variant_sparse_matrix$pos==13352&Donor2_variant_sparse_matrix$vaf>0,]



grep "13352" GTACTAATCGCTCCAT-1.snv


# 一行代码提取非C>A/C>T突变的数字位置
positions <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                            mutations[!grepl("_C_[AT]$", mutations)]))
tmp<- Donor2_variant_sparse_matrix[Donor2_variant_sparse_matrix$pos%in%positions&Donor2_variant_sparse_matrix$vaf>0.1,]



overlap_5 <- unique(mutation_1)
overlap_5<- overlap_5[-which(overlap_5 %in% c(jinxu_ppl,human_ppl_v5_minreads1,human_ppl_v5_minreads3))]

2366_G_A
15657_T_C
15967_C_T
16137_A_G

final[final$Position%in%c(2366,15657,15967,16137),]

> final[final$Position%in%c(2366,15657,15967,16137),]
      Chrom Position Ref Cons Reads1 Reads2 VarFreq Strands1 Strands2 Qual1
21834  chrM    15657   T    Y     46      7  13.21%        2        2    36
26333  chrM    15967   C    Y     14      6     30%        2        2    37
28482  chrM    16137   A    R     26      4  13.33%        2        2    37
31164  chrM     2366   G    R      5      5     50%        2        2    37
      Qual2 Pvalue MapQual1 MapQual2 Reads1Plus Reads1Minus Reads2Plus
21834    37   0.98        1        1         25          21          3
26333    37   0.98        1        1          7           7          4
28482    37   0.98        1        1         12          14          2
31164    37   0.98        1        1          3           2          3
      Reads2Minus VarAllele                barcode
21834           4         C ATGACCAGTAACCACA-1.snv
26333           2         T CCGGTAGGTCGTAAAT-1.snv
28482           2         G TTCCCGCCAGGACCAA-1.snv
31164           2         A TGCACTTGTTACTTCA-1.snv

 less -S ./*/TGCACTTGTTACTTCA-1.counts|grep "2366"

2366%in% Donor1_variant_sparse_matrix$pos
2366%in% Donor2_variant_sparse_matrix$pos
2366%in% Donor3_variant_sparse_matrix$pos
2366%in% Donor4_variant_sparse_matrix$pos


less -S /md01/nieyg/project/mito_mutation/01_pipeline/07_ppl_v3/masked_SNVcalling_percell/TGCACTTGTTACTTCA-1.counts|grep "2366"


positions <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                            overlap_5))


tmp<- final[final$Position%in%positions,]

less -S ./*/CACTAAGGTTCGGTAA-1.snv|grep "11747"
less -S ./*/AGTTTGATCATCACTT-1.snv|grep "11747"
less -S ./*/GCCTATTGTAACCTAG-1.snv|grep "11915"
less -S ./*/TGCTTGTGTATCTGGA-1.snv|grep "11924"
less -S ./*/AAGACATAGGCGAAAC-1.snv|grep "12104"
less -S ./*/GCCACTAAGGAACCGG-1.snv|grep "12366"
less -S ./*/CTCTAAGCACTAGCGT-1.snv|grep "12476"
less -S ./*/TTGGCGGGTGCATTTC-1.snv|grep "12506"
less -S ./*/ACGTCCTTCCGGTATG-1.snv|grep "12523"
less -S ./*/AGGATATAGGACACTT-1.snv|grep "12559"



3847    11747 CACTAAGGTTCGGTAA-1.snv
3848    11747 AGTTTGATCATCACTT-1.snv
3886    11915 GCCTATTGTAACCTAG-1.snv
3887    11924 TGCTTGTGTATCTGGA-1.snv
3953    12104 AAGACATAGGCGAAAC-1.snv
5117    12366 GCCACTAAGGAACCGG-1.snv
6292    12476 CTCTAAGCACTAGCGT-1.snv
6298    12506 TTGGCGGGTGCATTTC-1.snv
6303    12523 ACGTCCTTCCGGTATG-1.snv
6311    12559 AGGATATAGGACACTT-1.snv
6316    12586 CTTCACTCAGGAACAT-1.snv
8075    12814 TAGCTTAAGCGGTTAT-1.snv
8076    12814 GTGATCAGTTAATGCG-1.snv


[nieyg@master output_minreads1]$  less -S ./*/TGCACTTGTTACTTCA-1.snv|grep "2366"
Chrom	Position	Ref	Cons	Reads1	Reads2	VarFreq	Strands1	Strands2	Qual1	Qual2	PvalueMapQual1	MapQual2	Reads1Plus	Reads1Minus	Reads2Plus	Reads2Minus	VarAllele
chrM	2366	G	R	5	5	50%	2	2	37	37	0.98	1	1	3	2	2	A


[nieyg@master output_minreads1]$ less -S ./*/AGGATATAGGACACTT-1.snv|grep "12559"
chrM	12559	C	Y	20	4	16.67%	2	2	37	37	0.98	1	1	11	92	2	T


[nieyg@master output_minreads1]$ head -n 2 ./*/TGCACTTGTTACTTCA-1.snv
Chrom	Position	Ref	Cons	Reads1	Reads2	VarFreq	Strands1	Strands2	Qual1	Qual2	PvalueMapQual1	MapQual2	Reads1Plus	Reads1Minus	Reads2Plus	Reads2Minus	VarAllele

only_jinxu_ppl <- unique(mutation_2)
only_jinxu_ppl<- only_jinxu_ppl[-which(only_jinxu_ppl %in% c(human_ppl_v1,human_ppl_v5_minreads1,human_ppl_v5_minreads3))]
positions <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                            only_jinxu_ppl))

somatic_snv[somatic_snv$position%in%positions,]

grep "388" /md01/jinxu/Project/mgatk-speedup/44_GenoByCell/variant_sparse_matrix.tsv | \
awk -F'\t' '$NF > 0' | awk -F'\t' '$7 > 4' |\
less -S

less -S ./*/ATTTGTGAGTTCCCAC-1.snv|grep "388"
ATTTGTGAGTTCCCAC-1

overlap_6 <- intersect(mutation_1,mutation_4)
overlap_6<- overlap_6[-which(overlap_6 %in% c(jinxu_ppl,human_ppl_v5_minreads1))]

8279_T_C

somatic_snv_v5

tmp3_1<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor1/Donor1_somatic.tsv",sep="\t")
tmp3_2<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor2/Donor2_somatic.tsv",sep="\t")
tmp3_3<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor3/Donor3_somatic.tsv",sep="\t")
tmp3_4<- read.csv("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor4/Donor4_somatic.tsv",sep="\t")

somatic_snv_v5_3<- rbind(tmp3_1,tmp3_2,tmp3_3,tmp3_4)

somatic_snv_v5_3[somatic_snv_v5_3$Position==8279,]

TCCATAAAGACAGGCG-1


grep "3707" Donor2_variant_sparse_matrix_withoutstrand.tsv | \
awk -F'\t' '$NF > 0' | awk -F'\t' '$7 > 4' |\
less -S



[nieyg@master PBMC_lib5_output]$ python /md01/nieyg/pipeline/mito_mutation/01_human/SNV_summary_v2.py -i /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/output/Donor2/ -o ./snv_test_minread3_donor2.tsv -v

开始处理所有细胞的SNV结果...
输入目录: /md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/output/Donor2/
输出文件: ./snv_test_minread3_donor2.tsv
最小深度: 2
文件模式: *.snv
生成压缩文件: 是
找到 4640 个snv文件
步骤1: 解析所有snv文件...
已解析 1000/4640 个文件
已解析 2000/4640 个文件
已解析 3000/4640 个文件
已解析 4000/4640 个文件
总变异记录数: 117030
总细胞数: 4548
步骤2: 确定每个位置的alt_base...
确定的突变数量: 1386
步骤3: 计算突变统计量...
已处理 1000/1386 个突变
保存结果到: ./snv_test_minread3_donor2.tsv
保存压缩版本到: ./snv_test_minread3_donor2.tsv.gz

=== 统计报告 ===
总突变数: 1166
平均VAF: 0.490555
平均confident细胞比例: 0.0217
平均VAF阳性细胞比例: 0.0220


