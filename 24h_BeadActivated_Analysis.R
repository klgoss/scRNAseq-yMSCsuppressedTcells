library(Seurat)
library(ggplot2)
library(dplyr)
library(escape)
library(dittoSeq)
library(ggpubr)
library(RCurl)
library(AnnotationHub)
library(rstatix)
library(viridis)
library(RColorBrewer)
library(msigdbr)
library(annotate)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)

# this script analyzes the 24hour transwell data - activated and suppressed. Gene expression data only.

set.seed(123)

# read in data and create Seurat object for activated
activated.data <- Read10X("activated/filtered_feature_bc_matrix/")
activated <- CreateSeuratObject(activated.data, project="activated", min.cells=3, min.features = 200)

# read in data and create Seurat object for suppressed
suppressed.data <- Read10X("suppressed/filtered_feature_bc_matrix/")
suppressed <- CreateSeuratObject(suppressed.data, project="suppressed", min.cells=3, min.features = 200)

merged <- merge(activated, y=suppressed, add.cell.ids = c("activated", "suppressed"))

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident")

merged <- subset(merged, subset = nCount_RNA < 20000 & nFeature_RNA < 4000 & percent.mt < 15)
VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident")

merged[["RNA"]] <- split(merged[["RNA"]], f=merged$orig.ident)

merged <- NormalizeData(merged) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(merged, ndims = 30)

merged <- FindNeighbors(merged, dims=1:20, reduction="pca")
merged <- FindClusters(merged, resolution=0.1, cluster.name="unintegrated_clusters")
merged <- RunUMAP(merged, dims=1:20, reduction="pca",reduction.name="umap.unintegrated")

DimPlot(merged, group.by = "orig.ident")
DimPlot(merged, label=T)

# cluster 3 isn't T cells - remove
merged <- subset(merged, idents='3', invert=T)

FeaturePlot(merged, features = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"), order=T)

merged <- JoinLayers(merged)
cluster_markers <- FindAllMarkers(merged, logfc.threshold = 1, only.pos=T, min.pct=0.1)
cluster_markers <- cluster_markers[cluster_markers$p_val_adj < 0.05,]

top10 <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

DoHeatmap(merged, features = top10$gene)

# NK cell markers
FeaturePlot(merged, features = c("FCGR3A", "NCAM1", "KLRD1", "KLRK1", "KLRC1", "NCR1"), order=T)

# cluster 2 is NK cells - remove
merged <- subset(merged, idents='2', invert=T)

# since we removed 2 clusters, renormalize/scale/run dim red
merged <- NormalizeData(merged) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(merged, ndims = 30)

merged <- FindNeighbors(merged, dims=1:20, reduction="pca")
merged <- FindClusters(merged, resolution=0.1, cluster.name="unintegrated_clusters")
merged <- RunUMAP(merged, dims=1:20, reduction="pca",reduction.name="umap.unintegrated")

#######################################################################
merged <- readRDS("2025_24hGEX.rds")

umap <- DimPlot(merged, group.by = "orig.ident", cols = c('#ED1E2D', '#00AEEF'))
ggsave(umap, filename = 'new_figures/24h_umap.svg')

merged$orig.ident <- ifelse(merged$orig.ident == "activated", "Control", "Suppressed")

Idents(merged) <- merged$orig.ident
#levels(merged) <-c("Suppressed", "Control")
#Idents(merged) <- factor(x=Idents(merged), levels = levels(merged))

merged <- JoinLayers(merged)
deg <- FindMarkers(merged, ident.1="Control", ident.2="Suppressed", group.by = "orig.ident", logfc.threshold = 0.25)
deg <- deg[deg$p_val_adj < 0.05,]
write.table(deg, file = "24h_DEG.tsv", quote=F, sep = "\t")

keyvals <- ifelse(
  deg$avg_log2FC < 0, '#00AEEF',
  ifelse(deg$avg_log2FC > 0, '#ED1E2D',
         'black'))

names(keyvals)[keyvals == '#00AEEF'] <- 'Suppressed'
names(keyvals)[keyvals == '#ED1E2D'] <- 'Control'

deg_volcano <- EnhancedVolcano::EnhancedVolcano(deg, pointSize=5, lab = rownames(deg), labSize=0, x='avg_log2FC',y = 'p_val_adj',FCcutoff = 0.25,pCutoff = 0.05, title = '24h control vs suppressed', subtitle = "", legendPosition = 'bottom', colCustom = keyvals, gridlines.major = F, gridlines.minor = F)
ggsave(deg_volcano, filename="24h_DEG_volcano.svg")

up_control <- deg[deg$avg_log2FC > 0,]
up_suppressed <- deg[deg$avg_log2FC < 0,]

control_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                        keys = rownames(up_control),
                                        columns = "ENTREZID",
                                        keytype = "GENENAME")

suppressed_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys = rownames(up_suppressed),
                                           columns = "ENTREZID",
                                           keytype = "GENENAME")


# Compare pathways
list <- list(control_entrez$ENTREZID, suppressed_entrez$ENTREZID)
names(list) <- c("Control", "Suppressed")

compGO <- compareCluster(geneCluster = list,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "bonferroni", readable = T)
go_dot <- dotplot(compGO, title = "GO Enrichment Analysis", by = "count") + scale_fill_viridis(option='plasma')


isr <- scan("FINAL_ISR.txt", what = "character")
isr_heat <- DoHeatmap(merged, features = isr, group.by = "orig.ident", group.colors = c('#ED1E2D', "#00AEEF")) + theme(axis.text.y=element_blank()) 


detox_genes <- c("GSTP1","TXN", "TXNRD1", "PRDX5", "GPX4", "RDH11", "CCS", "UCP2", "PLK3", "AIF1", "G6PD", "HSPA1A")
detox_dot <- DotPlot(merged, features = detox_genes) + RotatedAxis() + scale_color_viridis() + ylab("") + xlab("Gene") + coord_flip()


nrf2 <- VlnPlot(merged, features = "NFE2L2", group.by = "orig.ident", pt.size=0, cols = c('#ED1E2D', "#00AEEF")) + 
  geom_boxplot(width=0.3) +
  stat_compare_means(method="wilcox.test", label.x=1.45, label="p.format") +
  theme(legend.position="none")


################ ssGSEA with escape ###################
library(escape)

GS.hallmark <- getGeneSets(library="H", species = "Homo sapiens")

merged <- runEscape(merged, 
                    method = "ssGSEA",
                    gene.sets = GS.hallmark, 
                    groups = 1000, 
                    min.size = 5,
                    new.assay.name = "escape.H")


oxphos <- VlnPlot(merged, features = "HALLMARK-OXIDATIVE-PHOSPHORYLATION", assay="escape.H", pt.size=0, group.by = "orig.ident", cols = c('#ED1E2D', "#00AEEF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method="wilcox.test", label="p.format", label.x=1.3) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score')
ggsave(oxphos, filename="new_figures/24h_OXPHOS_Vln.svg")


glycolysis <- VlnPlot(merged, features = "HALLMARK-GLYCOLYSIS", assay="escape.H", pt.size=0, group.by = "orig.ident", cols = c('#ED1E2D', "#00AEEF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method="wilcox.test", label="p.format", label.x=1.35) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score')


GS.Kegg <- getGeneSets(library="C2", subcategory = "CP:KEGG", species = "Homo sapiens") 

merged <- runEscape(merged, 
                    method = "ssGSEA",
                    gene.sets = GS.Kegg, 
                    groups = 1000, 
                    min.size = 5,
                    new.assay.name = "escape.KEGG")

ppp <- VlnPlot(merged, features = "KEGG-PENTOSE-PHOSPHATE-PATHWAY", assay="escape.KEGG", pt.size=0, group.by = "orig.ident", cols = c('#ED1E2D', "#00AEEF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method="wilcox.test", label="p.format", label.x=1.35) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score')
ggsave(ppp, filename="new_figures/24h_PPP_Vln.svg")


