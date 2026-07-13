library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(dittoSeq)
library(escape)
library(scater)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

set.seed(123)

# creation of seurat object
# Load Control Data #
TG0 <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/actTG0/filtered_feature_bc_matrix/")
TG0_obj <- CreateSeuratObject(TG0, project = "TGO", min.cells = 3, min.features = 200)
TG0_obj[["group"]] <- "control"
TG0_obj[["phase"]] <- "G0"

TG1 <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/actTG1/filtered_feature_bc_matrix/")
TG1_obj <- CreateSeuratObject(TG1, project = "TG1", min.cells = 3, min.features = 200)
TG1_obj[["group"]] <- "control"
TG1_obj[["phase"]] <- "G1"

TSG2 <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/actTSG2/filtered_feature_bc_matrix/")
TSG2_obj <- CreateSeuratObject(TSG2, project = "TSG2", min.cells = 3, min.features = 200)
TSG2_obj[["group"]] <- "control"
TSG2_obj[["phase"]] <- "SG2"

# Load Experimental Data #
G0plusyMSC <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/G0plusyMSC/filtered_feature_bc_matrix/")
G0plusyMSC_obj <- CreateSeuratObject(G0plusyMSC, project = "G0plusyMSC", min.cells = 3, min.features = 200)
G0plusyMSC_obj[["group"]] <- "experimental"
G0plusyMSC_obj[["phase"]] <- "G0"

G1plusyMSC <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/G1plusyMSC/filtered_feature_bc_matrix/")
G1plusyMSC_obj <- CreateSeuratObject(G1plusyMSC, project = "G1plusyMSC", min.cells = 3, min.features = 200)
G1plusyMSC_obj[["group"]] <- "experimental"
G1plusyMSC_obj[["phase"]] <- "G1"

SG2plusyMSC <- Read10X(data.dir = "/Volumes/MyPassport/rotation1_data/SG2plusyMSC/filtered_feature_bc_matrix/")
SG2plusyMSC_obj <- CreateSeuratObject(SG2plusyMSC, project = "SG2plusyMSC", min.cells = 3, min.features = 200)
SG2plusyMSC_obj[["group"]] <- "experimental"
SG2plusyMSC_obj[["phase"]] <- "SG2"

# Combine datasets #
obj <- merge(TG0_obj, y = c(TG1_obj, TSG2_obj,G0plusyMSC_obj,G1plusyMSC_obj,SG2plusyMSC_obj), project = "combine")

# Features (for QC) #
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

VlnPlot(obj, features = "nFeature_RNA", pt.size = 0.3)
VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 2)

obj <- subset(obj, subset = nFeature_RNA < 7500 & nCount_RNA < 100000 & percent.mt < 25)

# Normalize (log norm- cpm) #
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e6)

# Find Variable Features #
obj <- FindVariableFeatures(obj, selection.method = "vst")
top10 <- head(VariableFeatures(obj), 10)

all.genes <- rownames(obj)

# Scale data and run PCA #
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

pca <- DimPlot(obj, reduction = "pca")

ElbowPlot(obj, ndims = 30) # use 20 dims

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20)

saveRDS(obj, "96h_BeadActivated_SeuratObj.rds")

# load Seurat object
obj <- readRDS("96h_BeadActivated_SeuratObj.rds")
umap <- DimPlot(obj, group.by = "orig.ident") + ggtitle("")

# Differential expression
degs <- FindMarkers(obj, ident.1 = "control", ident.2 = "experimental", group.by = "group", logfc.threshold = 0.25, min.pct=0.1)
degs <- degs[degs$p_val_adj < 0.05,]

# Gene ontology analysis
up_control <- degs[degs$avg_log2FC > 0,]
up_suppressed <- degs[degs$avg_log2FC < 0,]

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
df <- compGO@compareClusterResult 

oxidative_stress <- df[df$Description == "response to oxidative stress", ]$geneID 
oxidative_stress <- unlist(strsplit(oxidative_stress, "/"))

oxidative_stress_dot <- DotPlot(obj, features = oxidative_stress, group.by = "group", cols = c("purple", "orange")) + RotatedAxis() 
ggsave(oxidative_stress_dot, filename = 'oxidative_stress_dotplot.svg')

VlnPlot(obj, features = "NFE2L2", group.by = "orig.ident")

obj <- ScaleData(obj, features = oxidative_stress)
oxidative_stress_heat <- DoHeatmap(obj, features = oxidative_stress, group.by = "group") + theme_void()
ggsave(oxidative_stress_heat, filename = "oxidative_stress_heat.tiff")

cat(oxidative_stress, sep = "\n", file = "oxidative_stressgenes_forHeat.txt")

nrf2 <- VlnPlot(obj, features = "NFE2L2", group.by = "group") + theme(legend.position='none')
ggsave(nrf2, filename = "NFE2L2_vln.svg")

glycolysis <- df[df$ID=="GO:0061621",]$geneID
glycolysis_genes <- unlist(strsplit(glycolysis, "/"))

oxphos <- df[df$ID == "GO:0006119",]$geneID
oxphos_genes <- unlist(strsplit(oxphos, "/"))

metabolism_dot <- DotPlot(obj, features = c(glycolysis_genes, oxphos_genes), group.by = "group", cols = c('purple', 'orange')) + RotatedAxis()
ggsave("metabolism_dotplot.svg", metabolism_dot)




# ISR genes
ISR_genes <- scan("FINAL_ISR.txt", sep = ",", what = "character")

Idents(obj) <- obj$orig.ident
levels(obj) <- c("G0plusyMSC", "G1plusyMSC", "SG2plusyMSC", "TGO", "TG1", "TSG2")
a <- DoHeatmap(obj, features = ISR_genes, label = T, size = 3, group.by = "orig.ident") +
  theme(axis.text.y = element_blank())
