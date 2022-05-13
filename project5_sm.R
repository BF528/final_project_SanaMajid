##BF528 Project 5##
##Student: Sana Majid##


##Load in packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#tximport needs to be accessed through bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
library(tximport) 

##Project 4 Programmer role

#Set working directory
setwd("/projectnb/bf528/project_4_scrnaseq")

#load the salmon avelin counts file (GSM2230760_salmon_counts) into R
salmon <- file.path("GSM2230760__salmon_quant/alevin/quants_mat.gz")  
txi <- tximport(salmon, type="alevin")

#choose criteria to identify 
class(txi)
head(txi) #has $abundance, $counts, and $countsfromabundance (we need counts)
counts <- txi$counts

# Creating Seurat Object
panca <- CreateSeuratObject(counts = counts , min.cells = 3, min.features = 200)
panca 
#An object of class Seurat 
#26442 features across 15147 samples within 1 assay 
#Active assay: RNA (26442 features, 0 variable features)

# Data QC - Create Mitochondrial Percent
panca[["percent.mt"]] <- PercentageFeatureSet(panca, pattern = "^MT-")

# switch to directory where you want to load files and graphs
setwd("/projectnb/bf528/users/dachsund/project_5_sm/")

#generate the violin plot for QC metrics
vln_qc <- VlnPlot(panca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
vln_qc #saved img as vln_qc.png

#generate feature scatter plots that help visualize feature-feature relationships
plot1 <- FeatureScatter(panca, feature1 = "nCount_RNA", feature2 = "percent.mt")
#warning msg: In cor(x = data[, 1], y = data[, 2]) : the standard deviation is zero
plot2 <- FeatureScatter(panca, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 #saved as scatter_plots.png

panca_subset <- subset(panca, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)

#Normalizing the data
panca <- NormalizeData(panca, normalization.method = "LogNormalize", scale.factor = 10000)
#Performing log-normalization
#0%   10   20   30   40   50   60   70   80   90   100%
#[----|----|----|----|----|----|----|----|----|----|
#**************************************************|
#Identification of highly variable features
panca <- FindVariableFeatures(panca, selection.method = "vst", nfeatures = 2000)
#Calculating gene variances
#0%   10   20   30   40   50   60   70   80   90   100%
#[----|----|----|----|----|----|----|----|----|----|
#**************************************************|
#Calculating feature variances of standardized and clipped values
#0%   10   20   30   40   50   60   70   80   90   100%
#[----|----|----|----|----|----|----|----|----|----|
#**************************************************|

top10 <- head(VariableFeaturePlot(panca), 10)
top10
# plot variable features with and without labelsVizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
plot1 <- VariableFeaturePlot(pbmc)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE,xnudge =0,ynudge= 0)
plot1 + plot2 ##saved as var_featues_plot
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

#scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#run linear dimension reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #saved as vizdimload.png
DimPlot(pbmc, reduction = "pca") #saved as dimplot_1

#Determine the ‘dimensionality’ of the dataset
panca_jack <- JackStraw(pbmc, num.replicate = 100)
panca_js <- ScoreJackStraw(panca_jack, dims = 1:20)
JackStrawPlot(panca_js, dims = 1:15)

#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15) #saved as jackstrawplot.png



pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#Number of nodes: 15147
#Number of edges: 461646

#Running Louvain algorithm...
#0%   10   20   30   40   50   60   70   80   90   100%
#[----|----|----|----|----|----|----|----|----|----|
#**************************************************|
#Maximum modularity in 10 random starts: 0.8718
#Number of communities: 13
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap") #saved as dimplot
ElbowPlot(pbmc) #saved as elbowplot

# Save RDS
saveRDS(cells_subset, file = "output.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers)


# Cells Per Cluster
cluster_numbers <- table(pbmc@meta.data$seurat_clusters)
cluster_numbers <- as.data.frame(cluster_numbers)
cluster_numbers %>% ggplot(mapping = aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity") +
  ggtitle("Cells Across Clusters") +
  xlab("Cluster ID") +
  ylab("Number of Cells") #saved as cell_dist


##Project 4 Analyst role

setwd("/projectnb/bf528/project_4_scrnaseq")
df <- readRDS("GSM2230760_seurat.rda")
setwd("/projectnb/bf528/users/dachsund/project_5_sm/")
# Identify marker genes for each cluster
cluster_markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25 )
signif <- cluster_markers[cluster_markers$p_val_adj<0.05,]
write.csv(signif, file="cluster_markers.csv")

cluster0 <- FindMarkers(df, ident.1=0, min.pct = 0.25)
head(cluster0, n=5) #SST, MT-ND1, MT-ND4, EEF1A1
cluster1 <- FindMarkers(df, ident.1=1, min.pct = 0.25)
head(cluster1, n=5) #INS, GNAS, DLK1, GCG, CDKN1C
cluster2 <- FindMarkers(df, ident.1=2, min.pct = 0.25)
head(cluster2, n=5) #FN1, COL1A1, IGFBP5, MT-RNR2, GCG
cluster3 <- FindMarkers(df, ident.1=3, min.pct = 0.25)
head(cluster3, n=5) #CELA3A, REG1B, REG1A, SPINK1, CPA1
cluster4 <- FindMarkers(df, ident.1=4, min.pct = 0.25)
head(cluster4, n=5) #GCG, TTR, TM4SF4, CHGB, ALDH1A1
cluster5 <- FindMarkers(df, ident.1=5, min.pct = 0.25)
head(cluster5, n=5) #KRT19, MMP7, PMEPA1, ANXA2, TACSTD2
cluster6 <- FindMarkers(df, ident.1=6, min.pct = 0.25)
head(cluster6, n=5) #MTND2P28, SAMD11, SDF4, MRPL20, ATAD3C
cluster7 <- FindMarkers(df, ident.1=7, min.pct = 0.25)
head(cluster7, n=5) #SPATS1, ACER3, GPR82, RPS6KA5, HELLPAR
cluster8 <- FindMarkers(df, ident.1=8, min.pct = 0.25)
head(cluster8, n=5) #CFAP74, NBL1, WNT4, STMN1, SSX2IP
cluster9 <- FindMarkers(df, ident.1=9, min.pct = 0.25)
head(cluster9, n=5) #FMOD, PXDN, COL3A1, COL5A2, FSTL1
cluster10 <- FindMarkers(df, ident.1=10, min.pct = 0.25)
head(cluster10, n=5) #HES4, AGRN, ERRFI1, DHRS3, TMEM51
cluster11 <- FindMarkers(df, ident.1=11, min.pct = 0.25)
head(cluster11, n=5) #CTRC, CELA2A, CELA3B, CELA3A, TCEA3
cluster12 <- FindMarkers(df, ident.1=12, min.pct = 0.25)
head(cluster12, n=5) #C1QA, C1QC, C1QB, LAPTM5, CD53

# Label clusters as a cell type based on marker genes
df_labels<-signif%>%mutate(celltype =
                             case_when(signif$cluster == "0" ~ "Delta",
                                       signif$cluster == "1" ~ "Beta_1",
                                       signif$cluster == "2" ~ "Gamma",
                                       signif$cluster == "3" ~ "Acinar",
                                       signif$cluster == "4" ~ "Alpha",
                                       signif$cluster == "5" ~ "Ductal",
                                       signif$cluster == "6" ~ "Beta_2",
                                       signif$cluster == "7" ~ "Exocrine glandular_1",
                                       signif$cluster== "8" ~ "Epsilon cells",
                                       signif$cluster == "9" ~ "Endothelial",
                                       signif$cluster == "10" ~ "Macrophage_1",
                                       signif$cluster == "11" ~ "Exocrine Glandular_2",
                                       signif$cluster == "12" ~ "Macrophages_2"))

# Visualize the clustered cells using a projection method
df <- RenameIdents(object = df, '0' = "Delta",'1' = "Beta_1",'2' = "Gamma",'3' = "Acinar",'4' = "Alpha",'5' = "Ductal",'6' = "Beta_2",'7' = "Exocrine glandular_1",'8' = "Epsilon cells",'9' = 'Endothelial','10' = 'Macrophage_1','11' = 'Exocrine Glandular_2','12' = 'Macrophages_2')
DimPlot(df, reduction = "umap", label = 'true')
          

# Visualize the top marker genes per cluster
top5 <- signif %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pl <- DoHeatmap(df, features = top10$gene) 
pl

# novel markers: top 2 of each cluster
FeaturePlot(df, features = c("PRRG3", "DLK1","FN1","COL1A1","TTR","CXCL1","KRT19",'EEF1A2','EDN3','ACER3','PLCE1', 'SPARC','ALDOB','PRSS21'))
