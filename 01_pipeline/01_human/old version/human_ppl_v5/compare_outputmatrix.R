library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 步骤1: 读取原始变异数据
variant_file <- "/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor1/Donor1_variant_sparse_matrix_withoutstrand.tsv.gz"
variants <- fread(variant_file)
variants <- variants[which(variants$vaf>0),]
variants$snv <- paste0(variants$pos,paste(variants$ref_base,variants$alt_base,sep=">"))
#variants <- variants[which(variants$snv%in%mutation_list),]
variants_new<- variants
variant_file2 <- "/md01/jinxu/Project/mgatk-speedup/44_GenoByCell/variant_sparse_matrix.tsv"
variants <- fread(variant_file2)
variants <- variants[which(variants$vaf>0),]
variants$snv <- paste0(variants$pos,paste(variants$ref_base,variants$alt_base,sep=">"))
#variants <- variants[which(variants$snv%in%mutation_list),]
variants_old<- variants

# for a cell
variants_new<- as.data.frame(variants_new)
variants_old<- as.data.frame(variants_old)

cell="AGGACGTAGGCTGTGC-1"

snv_new<- fread("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/donor_results/Donor1/Donor1_snv.tsv.gz")
snv_old<- fread("/md01/jinxu/Project/mgatk-speedup/13_coverge_pv/snv.tsv")

cell_snv<- fread("/md01/nieyg/project/mito_mutation/03_human_PBMC/PBMC_lib5_output/output/Donor1/AGGACGTAGGCTGTGC-1.snv")
cell_snv$Freq <- cell_snv$Reads2/(cell_snv$Reads2+cell_snv$Reads1)*100
cell_snv[cell_snv$Freq<90,]
cell_snv[cell_snv$Position%in%c(11600,14643,7897)]

Chrom   Position        Ref     Cons    Reads1  Reads2  VarFreq Strands1        Strands2        Qual1   Qual2   Pvalue  M

chrM    11600   G       R       37      2       5.13%   2       2       36      37      0.98    1       1       18      1


snv_new[snv_new$Position%in%c(210,5460,13395)]



mutation_list<- c("5460G>A","11600G>A","14643C>T","7897G>A")

variants_new_cell<- variants_new[which(variants_new$cell%in%cell),]

variants_old_cell<- variants_old[which(variants_old$cell==cell),]



# 步骤2: 读取细胞类型信息
celltype_file <- "/data/R01/renwx5/MT_mution/mutation_outs/sample_annotation/Lib5_barcode_sample_annotation.csv"
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





