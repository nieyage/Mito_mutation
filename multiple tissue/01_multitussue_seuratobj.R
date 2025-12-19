
conda activate r_cellchat 

library(BSgenome)
library(Matrix)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(hdf5r)
library(EnsDb.Mmusculus.v79)

#library(DropletUtils)

set.seed(1234)

# create a Seurat object containing the filtered RNA and ATAC data;
# 1. RNA filtered matrix 
counts <- Read10X_h5("/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/atac_fragments_filtered.tsv.gz"
filtered_RNA_data<-read.table("/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/counts_no_double_umi_001.tsv.gz")
filtered_RNA_counts = spread(filtered_RNA_data, V3, V2)
filtered_RNA_counts[is.na(filtered_RNA_counts)] <- 0
rownames(filtered_RNA_counts)<-filtered_RNA_counts$V1
filtered_RNA_counts<-filtered_RNA_counts[,-1]
RNA_barcode<-paste(colnames(filtered_RNA_counts),"-1",sep="")
colnames(filtered_RNA_counts)<-RNA_barcode

# 2. ATAC filtered matrix 

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
region_positive <- read.table(
  file = "/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
gr <- makeGRangesFromDataFrame(region_positive)


filtered_frag <- CreateFragmentObject(
  path = fragpath,
  cells = RNA_barcode,
  validate.fragments = FALSE
)
filtered_chrom_matrix<-FeatureMatrix(
  fragments = filtered_frag,
  features = gr
)
filtered_chrom_assay <- CreateChromatinAssay(
  counts = filtered_chrom_matrix,
  sep = c(":", "-"),
  fragments = filtered_frag,
  annotation=annotation
)

positive <-  CreateSeuratObject(counts = filtered_RNA_counts, min.cells=3, assay = "RNA")
positive[["ATAC"]]<- filtered_chrom_assay
objList<- list(positive)
# calculate the score of NS and TSS
for (i in seq_len(length(objList))) {
  DefaultAssay(objList[[i]]) <- "ATAC"
    objList[[i]] <- NucleosomeSignal(objList[[i]])
    objList[[i]] <- TSSEnrichment(objList[[i]],fast=FALSE)
    }
library(DataVisualizations)

# pdf("./01_QC/DensityScatter_QC.pdf")
# DensityScatter(objList[[1]], x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
# dev.off()

pdf("./01_QC/TSS_distribution_Positive.pdf")
  objList[[1]]$high.tss<-ifelse(objList[[1]]$TSS.enrichment > 4, 'High', 'Low')
  TSS<-TSSPlot(objList[[1]], group.by = 'high.tss') + NoLegend()+ labs(title = "NE")
  objList[[1]]$nucleosome_group <- ifelse(objList[[1]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
   print(TSS);
dev.off();

pdf("./01_QC/FragmentHistogram_QC.pdf")
FragmentHistogram(objList[[1]], group.by = 'nucleosome_group', region = "chr1-1-20000000")
dev.off()

sample<- c("Positive")
objList[[1]][["percent.mt"]] <- PercentageFeatureSet(objList[[1]], pattern = "^mt-")


  pdf(file = paste("./01_QC/3b_QC_before_",sample[1],".pdf", sep = ""),width=10,height=6)
  ##qc<-VlnPlot(object = objList[[1]],
   #         features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","nFeature_RNA"),
   #         ncol = 5,
   #         pt.size = 0.01
   #       )
  densityRNA<-plot(density(objList[[1]]@meta.data$nCount_RNA),xlim=c(0,5000))
  densityRNA_feature<-plot(density(objList[[1]]@meta.data$nFeature_RNA),xlim=c(0,5000))
  densityATAC<-plot(density(objList[[1]]@meta.data$nCount_ATAC),xlim=c(0,1000))
  #print(qc)
  print(densityRNA)
  print(densityRNA_feature)
  print(densityATAC)
  dev.off();



   # filter out low quality cells
# To remove doublets,select different cutoff#####

positive_filter <- subset(objList[[1]],
     subset= nCount_ATAC < 2000 & nCount_ATAC > 15 &
      nCount_RNA < 12000 & nCount_RNA > 50 &
      nFeature_RNA < 10000 & nFeature_RNA > 50 &
      TSS.enrichment >= 1 &
      nucleosome_signal < 2 
      )


## Load the dataset
samples <- c("Positive")

## Quality control

  positive_filter[["percent.mt"]] = PercentageFeatureSet(positive_filter,pattern="^mt-")
  pdf("./01_QC/_qc_plot.pdf",width=9)
  #p1=VlnPlot(positive_filter,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=ggplot(data=positive_filter[["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
  p3=ggplot(data=positive_filter[["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
  p4=ggplot(data=positive_filter[["percent.mt"]],aes(x=percent.mt))+geom_density()
  p5=FeatureScatter(positive_filter, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6=FeatureScatter(positive_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()


objList<-objList2

# Peak Calling 
# quantify counts in each peak

    peaks <- CallPeaks(positive_filter, macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2")
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(positive_filter),
      features = peaks,
      cells = colnames(positive_filter)
    )

    # create a new assay using the MACS2 peak set and add it to the Seurat object
positive_filter[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = fragpath,
      annotation = annotation
    )

#######integrate RNA and ATAC#####################
# RNA analysis
DefaultAssay(positive_filter) <- "RNA"
positive_filter<- FindVariableFeatures(positive_filter)
positive_filter<-ScaleData(positive_filter)


OSN_mm <- RunPCA(positive_filter) 
pdf("./02_All_celltype/ElbowPlot_RNA.pdf")
ElbowPlot(OSN_mm,ndims=50,reduction='pca')
dev.off();

OSN_mm <- RunUMAP(OSN_mm, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(OSN_mm) <- "ATAC"
OSN_mm <- RunTFIDF(OSN_mm)
OSN_mm <- FindTopFeatures(OSN_mm, min.cutoff = 20)
OSN_mm <- RunSVD(object=OSN_mm)


pdf("./02_All_celltype/LSI-depth-correalation.pdf")
DepthCor(OSN_mm)
dev.off();


# build a joint neighbor graph using both assays
OSN_mm <- FindMultiModalNeighbors(
  object = OSN_mm,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:40),
  #modality.weight.name = "RNA.weight",
  verbose = TRUE
)
OSN_mm <- RunUMAP(OSN_mm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
OSN_mm <- FindClusters(OSN_mm, graph.name = "wsnn", resolution =3, algorithm = 3, verbose = FALSE)

###reorder the level of sample#####

my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./02_All_celltype/OSN_mm_cluster_WNN.pdf",width=8,height=8)
###cluster
DimPlot(OSN_mm, cols=my47colors,pt.size=0.5, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
dev.off()


###### marker gene annotation ######
# RNA tree 

Idents(OSN_mm)<- OSN_mm$seurat_clusters
object<- OSN_mm
embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

library(ape)
library(circlize)  # 需要安装：install.packages("circlize")

pdf("./02_All_celltype/WNN_cluster_RNA_tree.pdf", width =6, height = 6)

par(mar = rep(0, 4))
ape::plot.phylo(data.tree, 
                type = "fan",
                show.tip.label = TRUE,
                tip.color = my47colors,
                edge.width = 2,
                edge.color = "gray50",
                label.offset = 0.02,
                cex = 0.8,
                rotate.tree = 0)
dev.off()


# ATAC tree 

embeddings <- Embeddings(object = object, reduction = "lsi")[,2:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
pdf("./02_All_celltype/WNN_cluster_ATAC_tree.pdf", width =6, height = 6)

par(mar = rep(0, 4))
ape::plot.phylo(data.tree, 
                type = "fan",
                show.tip.label = TRUE,
                tip.color = my47colors,
                edge.width = 2,
                edge.color = "gray50",
                label.offset = 0.02,
                cex = 0.8,
                rotate.tree = 0)
dev.off()


#####Annotate cells by RNA assay################
DefaultAssay(OSN_mm) <- "RNA"
Idents(OSN_mm)<-OSN_mm$seurat_clusters
# visualize the activities of canonical marker genes
liver_markers <- list(
  "T cell" = c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1"),
  "Hepatoblast" = c("Afp", "Sox9", "Epcam", "Krt19"),
  "Progenitor cell" = c("Prom1", "Kit"),
  "Hepatocyte" = c("Alb", "Gys2", "Ttr", "Cyp3a11", "Apoa1", "Apoa2"),
  "B cell" = c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Cd22"),
  "Dendritic cell" = c("Itgax", "Cd80", "Cd86", "H2-Ab1", "H2-Aa"),
  "Endothelial cell" = c("Pecam1", "Cdh5", "Kdr", "Fcgr2b", "Lyve1", "Stab2"),
  "Hepatic stellate cell" = c("Acta2", "Des", "Gfap", "Lrat", "Pdgfrb"),
  "Kupffer cell" = c("Cd68", "Clec4f", "Adgre1", "Fcgr1", "Csf1r"),
  "Monocyte" = c("Itgam", "Ly6c2", "Ccr2", "Cx3cr1"),
  "Myofibroblast" = c("Acta2", "Des", "Pdgfra", "Pdgfrb"),
  "Natural killer cell" = c("Ncr1", "Klrb1c", "Klrc1", "Itga1", "Nkg7"),
"Macrophage" = c("C1qa", "C1qb", "Ms4a7", "Adgre1")
)

lung_markers <- list(
  "Macrophage" = c("C1qa", "C1qb", "Ms4a7", "Adgre1"), 
  "Myofibroblast" = c("Acta2", "Myh11", "Tagln"),
  "Clara cell (Club cell)" = c("Scgb1a1", "Scgb3a2", "Cyp2f2"),
  "Epithelial cell" = c("Cdh1", "Epcam"),
  "Lipofibroblast" = c("Plin2", "Fabp4", "Pparg"),
  "Airway Dendritic Cell" = c("Itgae", "Xcr1", "Batf3"),
  "Dendritic Cell" = c("Itgax", "Tlr9"),
  "Eosinophil" = c("Itgav", "Siglecf"),
  "Goblet cell" = c("Muc5ac"),
  "Mesothelial cell" = c("Msln", "Wt1", "Upk3b"),
  "Neuroendocrine cell" = c("Calca"), 
  "Neutrophil" = c("Ly6g"),
  "T cell" = c("Cd3d", "Cd4", "Thy1"),
  "Type I Pneumocyte" = c("Ager", "Pdpn"),
  "Type II Pneumocyte" = c("Sftpc", "Slc34a2", "Etv5", "Lamp3", "Abca3")
)

kidney_markers <- list(
"Endothelial cell (Endo)" = c("Kdr", "Pecam1", "Cdh5"),
"Podocyte (Podo)" = c("Nphs1", "Nphs2", "Cdkn1c", "Bcam"),
"Proximal Tubule (PT)" = c("Slc34a1", "Slc5a2", "Slc22a6", "Slc22a8", "Cubn"),
"Loop of Henle, Ascending (LOH)" = c("Slc12a1", "Umod", "Cldn16"),
"Distal Convoluted Tubule (DCT)" = c("Slc12a3", "Calb1", "Trpm6"),
"Collecting Duct Principal Cell (CD-PC)" = c("Aqp2", "Scnn1b", "Scnn1g", "Hsd11b2", "Nr3c2"),
"Collecting Duct Intercalated Cell (CD-IC)" = c("Atp6v1g3", "Atp6v1b1", "Slc4a1", "Foxi1"),
"Collecting Duct Transitional Cell (CD-Trans)" = c("Sy17", "Parm1", "Sec23b"),
"Fibroblast (Fib)" = c("Pdgfrb", "Col1a1", "Col3a1"),
"Macrophage (Macro)" = c("Adgre1", "C1qa", "C1qb", "Cd68"),
"Neutrophil (Neutro)" = c("Ly6g", "S100a8", "S100a9"),
"B lymphocyte (B lymph)" = c("Cd19", "Ms4a1", "Cd79a"),
"T lymphocyte (T lymph)" = c("Cd3d", "Cd3e", "Cd4", "Cd8a"),
"Natural Killer Cell (NK)" = c("Ncr1", "Klrb1c", "Nkg7"))

spleen_markers <- list(
"T cell" = c("Cd3e", "Cd3d", "Cd4", "Cd8a", "Cd8b"),
"NK cell" = c("Nkg7", "Ncr1", "Klrb1c", "Kirc1", "Gzma"),
"B cell" = c("Cd79a", "Ms4a1", "Cd19", "Cr2", "Ighm"),
"Red blood cell (RBC)" = c("Hba-a1", "Hbb-bs", "Hba-a2"),
"Macrophage (Macro)" = c("C1qb", "C1qa", "Cd68", "Adgre1", "Cd14"),
"Monocyte (Mono)" = c("Lyz2", "Itgam", "Ccr2", "Cx3cr1"),
"Myeloid dendritic cell (mDC)" = c("Ccr7", "Itgax", "Cd80", "Cd86"),
"Plasmacytoid dendritic cell (pDC)" = c("Siglech", "Bst2", "Irfr7"),
"Mast cell" = c("Cpa3", "Mcpt4", "Mcpt1"),
"Plasma cell" = c("Ichain", "Jchain", "Igkc", "Ighg"),
"Neutrophil (Neutro)" = c("S100a8", "S100a9", "Lyz2", "Cxcr2", "Ccl3")
)


brain_markers<- c("Slc17a7","Slc17a6",#Glu+neurons
             "Slc32a1","Gad1","Gad2",#GABA+neurons 
             "Gfap","Aqp4", "Aldoc", "Fgfr3", #astrocyte
             "Cx3cr1","Aif1","Tmem119","Csf1r",#microglia
             "Mbp","Plp1","Olig2","Mag",#Oligodendrocyte
             "Pdgfra","Cspg4", "Tnr" #opc
)

colon_markers <- c(
  "c-Kit","CD166","CD24","CD44","Lgr5", #Colonic stem cell
  "c-Myc","D44","CD133", #Stem cell
  "Cbln2","Fam19a1","Hpse","Mrgprd","Ntm","Smr2","Trpa1", #Colonic sensory neuron
  "Cbln2","Fam19a1","Hpse","Mrgprd","Ntm","Smr2","Trpa1", #Neuron
  "Fam19a1","Ldhb","NECAB2","Nefh","Spp1" #Colonic neuron
)

skin_markers <- c(
  "Adipoq","Camp","Cebpa","Cebpb","LEPR","MMP1","MMP3","PDGFRA","PREF1", #Myofibroblast
  "GR-1", #Neural stem cell
  "IL13","IL22", #Skin cell
  "Lgr6" #Skin progenitor
)
oe_markers <- c("Syt1", # neurons
              "Omp",# Mature OSNs
              "Nqo1","Ncam2", # dorsal/ventral OSNs
              "Gap43","Gng8",#Immature OSNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#Erythrocytes
              "Mcpt8","Ccl4", #Basophils
              "C1qa","Ms4a7"#Macrophages
)


##### Annotate cells by RNA assay ################
DefaultAssay(OSN_mm) <- "RNA"
Idents(OSN_mm) <- OSN_mm$seurat_clusters

# 合并所有marker列表（按组织来源）
all_markers <- c(
  liver_markers,
  lung_markers, 
  kidney_markers,
  spleen_markers,
  list(
    "Brain_cell" = brain_markers,
    "Colon_cell" = colon_markers,
    "Skin_cell" = skin_markers,
    "OE_cell" = oe_markers
  )
)

# 展开所有基因
all_genes <- unique(unlist(all_markers))

# 检查哪些基因在数据中存在
available_genes <- intersect(all_genes, rownames(OSN_mm))

cat("总共", length(all_genes), "个marker基因\n")
cat("在数据中找到", length(available_genes), "个基因\n")

if (length(available_genes) < length(all_genes)) {
  missing_genes <- setdiff(all_genes, rownames(OSN_mm))
  cat("缺失", length(missing_genes), "个基因:\n")
  print(head(missing_genes, 20))
}

# 方法1：创建点图（DotPlot） - 按细胞类型分组
# 先按组织来源分别可视化

# 1. 肝脏marker基因
liver_available <- intersect(unlist(liver_markers), rownames(OSN_mm))
if (length(liver_available) > 0) {
  pdf("./02_All_celltype/liver_markers_dotplot.pdf", width = 12, height = 8)
  p <- DotPlot(OSN_mm, features = liver_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Liver Marker Genes Expression")
  print(p)
  dev.off()
}

# 2. 肺marker基因
lung_available <- intersect(unlist(lung_markers), rownames(OSN_mm))
if (length(lung_available) > 0) {
  pdf("./02_All_celltype/lung_markers_dotplot.pdf", width = 12, height = 10)
  p <- DotPlot(OSN_mm, features = lung_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Lung Marker Genes Expression")
  print(p)
  dev.off()
}

# 3. 肾脏marker基因
kidney_available <- intersect(unlist(kidney_markers), rownames(OSN_mm))
if (length(kidney_available) > 0) {
  pdf("./02_All_celltype/kidney_markers_dotplot.pdf", width = 14, height = 10)
  p <- DotPlot(OSN_mm, features = kidney_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Kidney Marker Genes Expression")
  print(p)
  dev.off()
}

# 4. 脾脏marker基因
spleen_available <- intersect(unlist(spleen_markers), rownames(OSN_mm))
if (length(spleen_available) > 0) {
  pdf("./02_All_celltype/spleen_markers_dotplot.pdf", width = 12, height = 10)
  p <- DotPlot(OSN_mm, features = spleen_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Spleen Marker Genes Expression")
  print(p)
  dev.off()
}

# 5. 嗅上皮marker基因（假设你的数据是嗅上皮）
oe_available <- intersect(oe_markers, rownames(OSN_mm))
if (length(oe_available) > 0) {
  pdf("./02_All_celltype/oe_markers_dotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = oe_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Olfactory Epithelium Marker Genes Expression")
  print(p)
  dev.off()
}

brain_available <- intersect(brain_markers, rownames(OSN_mm))
if (length(brain_available) > 0) {
  pdf("./02_All_celltype/brain_markers_dotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = brain_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Brain Marker Genes Expression")
  print(p)
  dev.off()
}

colon_available <- intersect(colon_markers, rownames(OSN_mm))
if (length(colon_available) > 0) {
  pdf("./02_All_celltype/colon_markers_dotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = colon_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Colon Marker Genes Expression")
  print(p)
  dev.off()
}

skin_available <- intersect(skin_markers, rownames(OSN_mm))
if (length(skin_available) > 0) {
  pdf("./02_All_celltype/skin_markers_dotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = skin_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Skin Marker Genes Expression")
  print(p)
  dev.off()
}



# 基于marker基因表达进行细胞注释
annotate_cells_by_markers <- function(seurat_obj, marker_lists) {
  # 为每个细胞类型计算marker基因的平均表达
  avg_exp <- matrix(0, nrow = length(marker_lists), 
                    ncol = length(levels(Idents(seurat_obj))))
  rownames(avg_exp) <- names(marker_lists)
  colnames(avg_exp) <- levels(Idents(seurat_obj))
  
  for (cell_type in names(marker_lists)) {
    markers <- marker_lists[[cell_type]]
    # 只取在数据中存在的基因
    markers <- intersect(markers, rownames(seurat_obj))
    
    if (length(markers) > 0) {
      # 计算该细胞类型marker基因在每个cluster的平均表达
      exp_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
      for (cluster in colnames(avg_exp)) {
        cells <- WhichCells(seurat_obj, idents = cluster)
        if (length(cells) > 0) {
          cluster_exp <- exp_data[markers, cells, drop = FALSE]
          avg_exp[cell_type, cluster] <- mean(as.matrix(cluster_exp))
        }
      }
    }
  }
  
  # 找到每个cluster表达最高的细胞类型
  cluster_annotations <- apply(avg_exp, 2, function(x) {
    if (all(x == 0)) return("Unknown")
    names(which.max(x))
  })
  
  return(list(
    avg_expression = avg_exp,
    annotations = cluster_annotations
  ))
}

# 应用注释
annotation_result <- annotate_cells_by_markers(OSN_mm, all_markers)

# 查看注释结果
cat("Cluster annotations:\n")
print(annotation_result$annotations)

# 保存注释结果
write.csv(data.frame(
  Cluster = names(annotation_result$annotations),
  Annotation = annotation_result$annotations,
  stringsAsFactors = FALSE
), "./02_All_celltype/cluster_annotations.csv", row.names = FALSE)

# 将注释添加到seurat对象
OSN_mm$celltype_annotation <- annotation_result$annotations[as.character(OSN_mm$seurat_clusters)]

# 可视化注释结果
pdf("./02_All_celltype/celltype_annotation_umap.pdf", width = 10, height = 8)
DimPlot(OSN_mm, 
        group.by = "celltype_annotation", 
        label = TRUE,
        repel = TRUE,
        pt.size = 0.5,reduction = "wnn.umap") +
  ggtitle("Cell Type Annotations")
dev.off()


# ATAC
library(BSgenome.Mmusculus.UCSC.mm10)
# ORN recall peak for subcluster 
DefaultAssay(OSN_mm)<-"ATAC"
peak<-CallPeaks(
       OSN_mm,
       group.by = "seurat_clusters",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 1.87e+09, #mm10 size
       outdir="./02_All_celltype/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(OSN_mm),
     features = peak,
     cells = colnames(OSN_mm)
     )     
OSN_mm[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(OSN_mm),
  annotation = Annotation(OSN_mm)
)

DefaultAssay(OSN_mm) <- 'peaks'
gene.activities <- GeneActivity(OSN_mm)
OSN_mm[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
OSN_mm <- NormalizeData(
  object = OSN_mm,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(OSN_mm$nCount_ACTIVITY)
)

DefaultAssay(OSN_mm) <- 'ACTIVITY'


# 1. 肝脏marker基因
liver_available <- intersect(unlist(liver_markers), rownames(OSN_mm))
if (length(liver_available) > 0) {
  pdf("./02_All_celltype/liver_markers_ACTIVITYdotplot.pdf", width = 12, height = 8)
  p <- DotPlot(OSN_mm, features = liver_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Liver Marker Genes Expression")
  print(p)
  dev.off()
}

# 2. 肺marker基因
lung_available <- intersect(unlist(lung_markers), rownames(OSN_mm))
if (length(lung_available) > 0) {
  pdf("./02_All_celltype/lung_markers_ACTIVITYdotplot.pdf", width = 12, height = 10)
  p <- DotPlot(OSN_mm, features = lung_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Lung Marker Genes Expression")
  print(p)
  dev.off()
}

# 3. 肾脏marker基因
kidney_available <- intersect(unlist(kidney_markers), rownames(OSN_mm))
if (length(kidney_available) > 0) {
  pdf("./02_All_celltype/kidney_markers_ACTIVITYdotplot.pdf", width = 14, height = 10)
  p <- DotPlot(OSN_mm, features = kidney_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Kidney Marker Genes Expression")
  print(p)
  dev.off()
}

# 4. 脾脏marker基因
spleen_available <- intersect(unlist(spleen_markers), rownames(OSN_mm))
if (length(spleen_available) > 0) {
  pdf("./02_All_celltype/spleen_markers_ACTIVITYdotplot.pdf", width = 12, height = 10)
  p <- DotPlot(OSN_mm, features = spleen_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Spleen Marker Genes Expression")
  print(p)
  dev.off()
}

# 5. 嗅上皮marker基因（假设你的数据是嗅上皮）
oe_available <- intersect(oe_markers, rownames(OSN_mm))
if (length(oe_available) > 0) {
  pdf("./02_All_celltype/oe_markers_ACTIVITYdotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = oe_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Olfactory Epithelium Marker Genes Expression")
  print(p)
  dev.off()
}

brain_available <- intersect(brain_markers, rownames(OSN_mm))
if (length(brain_available) > 0) {
  pdf("./02_All_celltype/brain_markers_ACTIVITYdotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = brain_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Brain Marker Genes Expression")
  print(p)
  dev.off()
}

colon_available <- intersect(colon_markers, rownames(OSN_mm))
if (length(colon_available) > 0) {
  pdf("./02_All_celltype/colon_markers_ACTIVITYdotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = colon_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Colon Marker Genes Expression")
  print(p)
  dev.off()
}

skin_available <- intersect(skin_markers, rownames(OSN_mm))
if (length(skin_available) > 0) {
  pdf("./02_All_celltype/Skin_markers_ACTIVITYdotplot.pdf", width = 16, height = 12)
  p <- DotPlot(OSN_mm, features = skin_available, 
               cols = c("lightblue", "red"), 
               dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Skin Marker Genes Expression")
  print(p)
  dev.off()
}




saveRDS(OSN_mm,"./02_All_celltype/WNN_integrated_all_celltype_recallpeaks.rds")

# Trackplot for annotation 
##Track for Marker genes promoters
Idents(OSN_mm)<-OSN_mm$seurat_clusters
library(BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(OSN_mm) <- "peaks"
# first compute the GC content for each peak
# conda r4-base环境有问题，用 r43计算
#conda activate r43
#.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r43/lib/R/library")
library(Matrix)
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

OSN_mm <- RegionStats(OSN_mm, genome = BSgenome.Mmusculus.UCSC.mm10)

#Annotation(OSN_mm)$tx_id <- Annotation(honeybee)$gene_name 


######Visulize track and RNA exp######
idents.plot <- Idents(OSN_mm)

 all_markers_list <- list(
    liver_markers,
    lung_markers,
    kidney_markers,
    spleen_markers
  )
  
  # 添加向量格式的marker
  all_markers_list <- c(all_markers_list, list(
    "Brain" = brain_markers,
    "Colon" = colon_markers,
    "Skin" = skin_markers,
    "OE" = oe_markers
  ))
  
  # 展开所有marker
  combined_markers <- list()
  for (marker_set in all_markers_list) {
    if (is.list(marker_set)) {
      combined_markers <- c(combined_markers, marker_set)
    } else if (is.character(marker_set)) {
      # 处理向量格式的marker
      combined_markers <- c(combined_markers, list("Custom" = marker_set))
    }
  }
combined_markers<- unlist(combined_markers)  
# link peaks to genes
OSN_mm <- LinkPeaks(
  object = OSN_mm,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = combined_markers
)
DefaultAssay(OSN_mm) <-"RNA"
input<- unique(combined_markers)[which(unique(combined_markers)%in%rownames(OSN_mm))]
DefaultAssay(OSN_mm) <-"peaks"
pdf("./02_All_celltype/Marker_gene-peaktrack-RNAexp-WNN.pdf",height=8,width=8)
for(i in input[158:180]){
  print(i)
  p1 <- CoveragePlot(
  object = OSN_mm,
  region = i,
  #features = i,
  #expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  annotation=TRUE,
  extend.downstream = 500
)
print(p1)}
dev.off()


Idents(OSN_mm)<-OSN_mm$seurat_clusters
#####further annotation########
OSN_mm <- RenameIdents(
  object = OSN_mm,
  '0' = 'T cell',
  '1' = 'Natural killer cell(NK)',
  '2' = 'Natural killer cell(NK)',
  '3' = 'Fibroblast(FB)',
  '4' = 'Proximal Tubule(PT)',
  '5' = 'OSN',
  '6' = 'Natural killer cell(NK)',
  '7' = 'T cell',
  '8' = 'Macrophage(MP)',
  '9' = 'Lipofibroblast(LipoFB)',
  '10' = 'Proximal Tubule(PT)',
  '11' = 'OSN',
  '12' = 'Neutrophil(Neutro)',
  '13' = 'Collecting Duct Principal Cell(CD-PC)',
  '14' = 'Proximal Tubule(PT)',
  '15' = 'B cell',
  '16' = 'Macrophage(MP)',
  '17' = 'Neutrophil(Neutro)',
  '18' = 'Myofibroblast(MyoFB)',
  '19' = 'Red blood cell(RBC)',
  '20' = 'Distal Convoluted Tubule(DCT)',
  '21' = 'Dendritic cell(DC)',
  '22' = 'OSN',
  '23' = 'Fibroblast(FB)',
  '24' = 'Plasma cell',
  '25' = 'Endothelial cell(EC)',
  '26' = 'T cell',
  '27' = 'Collecting Duct Principal Cell(CD-PC)',
  '28' = 'Hepatocyte',
  '29' = 'Proximal Tubule(PT)',
  '30' = 'Clara cell(Club cell)',
  '31' = 'OSN',
  '32' = 'Macrophage(MP)',
  '33' = 'Endothelial cell(EC)',
  '34' = 'Colonic sensory neuron',
  '35' = 'Endothelial cell(EC)',
  '36' = 'Proximal Tubule(PT)',
  '37' = 'Hepatocyte',
  '38' = 'Loop of Henle, Ascending(LOH)'
  )
OSN_mm@meta.data$cell_type<-Idents(OSN_mm)
order<- c("T cell","Natural killer cell(NK)","B cell","Macrophage(MP)","Dendritic cell(DC)","Neutrophil(Neutro)",
  "Fibroblast(FB)","Myofibroblast(MyoFB)","Lipofibroblast(LipoFB)",
  "Endothelial cell(EC)",
  "OSN",
  "Hepatocyte",
  "Proximal Tubule(PT)","Loop of Henle, Ascending(LOH)","Distal Convoluted Tubule(DCT)","Collecting Duct Principal Cell(CD-PC)","Red blood cell(RBC)","Plasma cell",
  "Clara cell(Club cell)","Colonic sensory neuron"
  )
OSN_mm@meta.data$cell_type<- factor(OSN_mm@meta.data$cell_type,levels=order)

table(OSN_mm$cell_type,OSN_mm$orig.ident)

Idents(OSN_mm)<- OSN_mm$cell_type

# Dimplot with annotation
pdf("./02_All_celltype/Multitissue_annotation_allcelltype_UMAP.pdf",width=8,height=5)
DimPlot(OSN_mm, label = T, repel = TRUE,pt.size=0.5, cols=my47colors, reduction = "wnn.umap",group.by = "cell_type")
DimPlot(OSN_mm, label = F, repel = TRUE,pt.size=0.5, cols=my47colors, reduction = "wnn.umap",group.by = "cell_type")+ ggtitle("")
dev.off()

# Save rds have annotation information 
DefaultAssay(OSN_mm) <- "RNA"
saveRDS(OSN_mm,"./02_All_celltype/WNN_integrated_all_celltype_with_annotation.rds")

OSN_mm<- readRDS("./02_All_celltype/WNN_integrated_all_celltype_with_annotation.rds")
# Dotplot with annotation
DefaultAssay(OSN_mm) <- 'RNA'
marker_list <- list(
  "T cell" = c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1", "Thy1"),
  "Natural killer cell(NK)" = c("Ncr1", "Klrb1c", "Klrc1", "Nkg7", "Gzma"),
  "B cell" = c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Cd22", "Cr2"),
  "Macrophage(MP)" = c("C1qa", "C1qb", "Ms4a7", "Adgre1", "Cd68", "Cd14", "Csf1r"),
  "Dendritic cell(DC)" = c("Itgax", "Cd80", "Cd86", "H2-Ab1", "H2-Aa",  "Ccr7"),
  "Neutrophil(Neutro)" = c("S100a8", "S100a9"),
  "Fibroblast(FB)" = c( "Col3a1","Pdgfra","Pdgfrb","Col1a1"),
  "Myofibroblast(MyoFB)" = c("Myh11", "Tagln"),
  "Lipofibroblast(LipoFB)" = c("Plin2", "Pparg"),
  "Endothelial cell(EC)" = c("Stab2","Cdh5", "Kdr", "Fcgr2b"),
  "OSN" = c("Syt1", "Omp"),
  "Hepatocyte" = c("Ttr"),
  "Proximal Tubule(PT)" = c("Slc34a1", "Slc5a2", "Slc22a6", "Slc22a8", "Cubn"),
  "Loop of Henle, Ascending(LOH)" = c("Slc12a1", "Umod", "Cldn16"),
  "Distal Convoluted Tubule(DCT)" = c("Slc12a3",  "Trpm6"),
  "Collecting Duct Principal Cell(CD-PC)" = c("Nr3c2" ,"Hsd11b2"),
  "Red blood cell(RBC)" = c( "Hbb-bs", "Hba-a2", "Hbb-bt"),
  "Plasma cell" = c("Jchain", "Igkc"),
  "Clara cell(Club cell)" = c( "Cyp2f2"),
  "Colonic sensory neuron" = c("Ntm")
)
markers<- unlist(marker_list)
library(RColorBrewer)
pdf("./02_All_celltype/OSN_mm_celltype-RNA_annotation.pdf",width=15,height=6)
p<-DotPlot(OSN_mm,cols = c("lightblue", "red"),  features = unique(markers),dot.scale = 6)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
DotPlot_df <- p$data
ggplot(data=DotPlot_df,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),)+
  scale_size_continuous(range=c(1,6),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_color_gradientn(colours =c(brewer.pal(9,"Blues")[2],brewer.pal(9,"Reds")[c(2,4,6)]) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",color="Scaled average expression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),axis.title=element_blank(),legend.position="bottom",legend.direction="horizontal")
dev.off()

DefaultAssay(OSN_mm) <- 'ACTIVITY'

marker_list <- list(
  "T cell" = c("Cd4", "Cd3e", "Cd3g", "Cd8b1", "Thy1"),
  "Natural killer cell(NK)" = c("Ncr1", "Klrc1",  "Nkg7"),
  "B cell" = c("Cd79a", "Cd79b", "Cd22", "Cr2", "Ighm"),
  "Macrophage(MP)" = c("C1qb", "Ms4a7"),
  "Dendritic cell(DC)" = c("Itgax", "Cd86", "H2-Aa"),
  "Neutrophil(Neutro)" = c("S100a8", "S100a9"),
  "Fibroblast(FB)" = c("Pdgfrb", "Col1a1"),
  "Myofibroblast(MyoFB)" = c("Myh11", "Tagln"),
  "Lipofibroblast(LipoFB)" = c("Plin2", "Fabp4", "Pparg"),
  "Endothelial cell(EC)" = c("Pecam1","Stab2"),
  "OSN" = c("Syt1", "Omp"),
  "Hepatocyte" = c("Alb", "Gys2", "Ttr"),
  "Proximal Tubule(PT)" = c("Slc34a1", "Slc5a2", "Slc22a6", "Slc22a8"),
  "Loop of Henle, Ascending(LOH)" = c("Slc12a1", "Cldn16"),
  "Distal Convoluted Tubule(DCT)" = c("Slc12a3", "Calb1", "Trpm6"),
  "Collecting Duct Principal Cell(CD-PC)" = c("Nr3c2", "Scnn1g" ),
  "Red blood cell(RBC)" = c("Hbb-bt"),
  "Plasma cell" = c(),
  "Clara cell(Club cell)" = c( "Cyp2f2"),
  "Colonic sensory neuron" = c("Mrgprd", "Ntm")
)
markers<- unlist(marker_list)
markers<- markers[which(markers%in%rownames(OSN_mm))]
library(RColorBrewer)

pdf("./02_All_celltype/OSN_mm_celltype-ATAC_annotation.pdf",width=15,height=6)
p<-DotPlot(OSN_mm,cols = c("lightblue", "#08306B"),  features = unique(markers),dot.scale = 6)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
#p
DotPlot_df <- p$data
ggplot(data=DotPlot_df,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),)+
  scale_size_continuous(range=c(1,9),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_color_gradientn(colours =c(brewer.pal(9,"Blues")) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",color="Scaled average expression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),axis.title=element_blank(),legend.position="bottom",legend.direction="horizontal")
dev.off()

OSN_mm[["percent.mt_RNA"]] = PercentageFeatureSet(OSN_mm,pattern="^mt-")

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

bam_file <- "/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/atac_possorted_bam.dedup.bam"
output_file <- "barcode_mtDNA_stats_from_bam.csv"

# 运行
stats <- extract_bam_stats("/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/atac_possorted_bam.dedup.bam", "/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/barcodes.txt")
print(head(stats, 10))
fwrite(stats, "barcode_mtDNA_stats_from_bam.csv")

write.csv(OSN_mm@meta.data,"/md01/nieyg/project/multitissue_plog/10x/multitissue_10x_masked/outs/barcode_metadata.csv")

OSN_mm$mtDNA_content <- stats[match(colnames(OSN_mm),stats$barcode),]$mt_ratio


pdf("Celltype_QC.pdf",width=10,height=25)
VlnPlot(object = OSN_mm,
           features = c("nCount_RNA","nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","mtDNA_content"),
           ncol = 1,
           pt.size = 0.01
         )
dev.off()


# 保留每个组织中比例较大的主要细胞类型标记
tissue_markers <- list()

# 1. 肝脏 - 主要细胞类型
tissue_markers$Liver <- list(
  "Hepatocyte" = c("Alb", "Gys2", "Ttr", "Cyp3a11", "Apoa1", "Apoa2")  # 肝细胞（主要功能细胞）
)
# 2. 肺部 - 主要细胞类型
tissue_markers$Lung <- list(
  "Pneumocyte" = c("Sftpc", "Slc34a2", "Abca3","Ager", "Pdpn") # 肺泡细胞（主要）
)
# 3. 肾脏 - 主要细胞类型
tissue_markers$Kidney <- list(
  "Proximal Tubule (PT)" = c("Slc34a1", "Slc5a2", "Slc22a6", "Cubn"),   # 近端小管（主要）
  "Distal Convoluted Tubule (DCT)" = c("Slc12a3", "Calb1", "Trpm6")   # 远端小管                    
)
# 4. 脾脏 - 主要细胞类型
tissue_markers$Spleen <- list(
  "Plasma cell" = c("Ichain", "Jchain", "Igkc", "Ighg"),                     # T细胞
  "Red blood cell (RBC)" = c("Hba-a1", "Hbb-bs", "Hba-a2")            # 红细胞
)
# 5. 大脑 - 主要细胞类型
tissue_markers$Brain <- list(
  "Neurons" = c("Slc17a7", "Slc17a6", "Slc32a1", "Gad1", "Gad2")      # 神经元（主要）
)

# 6. 嗅上皮 - 主要细胞类型
tissue_markers$Olfactory_Epithelium <- list(
  "Mature OSNs" = c("Omp")                                          # 成熟嗅觉神经元（主要）
)

# 7. 结肠 - 主要细胞类型
tissue_markers$Colon <- list(
  "Colonic Epithelial cells" = c("Lgr5", "CD44")                              # 上皮细胞（主要）
)

# 8. 皮肤 - 主要细胞类型
tissue_markers$Skin <- list(
  "Skin cell" = c("IL13","IL22")                               # 成纤维细胞
)


# 首先提取所有标记基因并添加组织标签
all_markers <- data.frame()

# 收集所有标记基因
for(tissue_name in names(tissue_markers)) {
  for(cell_type in names(tissue_markers[[tissue_name]])) {
    genes <- tissue_markers[[tissue_name]][[cell_type]]
    tissue_markers_df <- data.frame(
      Gene = genes,
      Tissue = tissue_name,
      CellType = cell_type,
      Label = paste0(tissue_name, "\n", cell_type),
      stringsAsFactors = FALSE
    )
    all_markers <- rbind(all_markers, tissue_markers_df)
  }
}

# 去除重复的基因（如果有）
all_markers <- all_markers[!duplicated(all_markers$Gene), ]

# 创建基因标签映射
gene_labels <- setNames(all_markers$Label, all_markers$Gene)

# 方法1：使用DotPlot（推荐）
library(Seurat)
library(ggplot2)

# 提取所有基因
genes_to_plot <- all_markers$Gene

# 检查哪些基因在数据中存在
available_genes <- genes_to_plot[genes_to_plot %in% rownames(OSN_mm)]

if(length(available_genes) > 0) {
  # 创建点图
  p <- DotPlot(OSN_mm, 
               features = available_genes,
               group.by = "cell_type",
               cols = c("lightgrey", "red"),
               dot.scale = 6) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    ggtitle("Tissue-Specific Marker Genes Expression") +
    labs(x = "Genes (Colored by Tissue)", y = "Clusters")
  
  # 在x轴标签中添加组织信息
  # 创建颜色映射
  tissue_colors <- c(
    Liver = "#E41A1C",      # 红色
    Lung = "#377EB8",       # 蓝色
    Kidney = "#4DAF4A",     # 绿色
    Spleen = "#984EA3",     # 紫色
    Brain = "#FF7F00",      # 橙色
    Olfactory_Epithelium = "#FFFF33",  # 黄色
    Colon = "#A65628",      # 棕色
    Skin = "#F781BF"        # 粉色
  )
  
  # 为每个基因添加颜色
  gene_colors <- setNames(tissue_colors[all_markers$Tissue[match(available_genes, all_markers$Gene)]], 
                          available_genes)
  
  # 自定义x轴标签
  custom_labels <- paste0(available_genes, "\n(", 
                         all_markers$Tissue[match(available_genes, all_markers$Gene)], ")")
  
  # 应用自定义标签
  p <- p + scale_x_discrete(labels = custom_labels)
  
  # 打印图形
  pdf("Tissue-Specific Marker Genes Expression(cell_type).pdf",width=14,height=6)
  print(p)
  dev.off()
} else {
  warning("No marker genes found in the dataset!")
}
