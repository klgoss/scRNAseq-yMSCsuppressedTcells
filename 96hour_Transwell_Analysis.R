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
library(rstatix)

set.seed(123)

# read in data
data.activated_top <- Read10X("activated_top/filtered_feature_bc_matrix/")
data.activated_bottom <- Read10X("activated_bottom/filtered_feature_bc_matrix/")
data.suppressed_top <- Read10X("suppressed_top/filtered_feature_bc_matrix/")
data.suppressed_bottom <- Read10X("suppressed_bottom/filtered_feature_bc_matrix/")

# Create Seurat objects
activated.top <- CreateSeuratObject(data.activated_top, project="activated.top", min.cells=3, min.features = 200)
activated.bottom <- CreateSeuratObject(data.activated_bottom, project="activated.bottom", min.cells=3, min.features = 200)
suppressed.top <- CreateSeuratObject(data.suppressed_top, project="suppressed.top", min.cells=3, min.features=200)
suppressed.bottom <- CreateSeuratObject(data.suppressed_bottom, project="suppressed.bottom", min.cells=3, min.features = 200)

# merge
combined <- merge(activated.top,y = c(activated.bottom, suppressed.top, suppressed.bottom))

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident")

combined <- subset(combined, subset = nCount_RNA < 40000 & nFeature_RNA < 6000 & percent.mt < 15)
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident")

combined[["RNA"]] <- split(combined[["RNA"]], f=combined$orig.ident)

combined <- NormalizeData(combined) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(combined, ndims = 30)

combined <- FindNeighbors(combined, dims=1:25, reduction="pca")
combined <- FindClusters(combined, resolution=0.1, cluster.name="unintegrated_clusters")
combined <- RunUMAP(combined, dims=1:25, reduction="pca",reduction.name="umap.unintegrated")

DimPlot(combined, group.by = "orig.ident")
DimPlot(combined, label=T)

combined <- IntegrateLayers(
  object = combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca")

combined <- FindNeighbors(combined, reduction = "integrated.rpca", dims = 1:25)
combined <- FindClusters(combined, resolution = 0.1, cluster.name = "rpca_clusters")
combined <- RunUMAP(combined, reduction = "integrated.rpca", dims = 1:25, reduction.name = "umap.rpca")

DimPlot(combined, reduction="pca", group.by = "orig.ident")

DimPlot(combined, reduction="umap.rpca")
DimPlot(combined, reduction="umap.rpca", group.by = "orig.ident")

Idents(combined) <- combined$orig.ident

combined <- JoinLayers(combined)

Idents(combined)
levels(combined)

new_order <- c("activated.bottom", "suppressed.bottom", "activated.top", "suppressed.top")
combined@active.ident <- factor(combined@active.ident, levels = new_order)

#saveRDS(combined, file = "2025_TopBottom.rds")

combined <- readRDS("2025_TopBottom.rds")


pca <- DimPlot(combined, reduction="pca", cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF"), shuffle=T)
ggsave(pca, filename="new_figures/2025_TopBottom_pca.svg")

markers <- FindAllMarkers(combined, logfc.threshold = 0.25, only.pos=T, min.pct = 0.1)
markers <- markers[markers$p_val_adj < 0.05,]
write.table(markers, file = "TopBottom_ALLMarkers.tsv", quote=F, sep = "\t")


txn <- VlnPlot(combined, features = "TXN", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) +
  geom_boxplot(width=0.3) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none")
ggsave("new_figures/96h_TopBottom_TXN_Vln.svg")

txn$data %>% dunn_test(TXN ~ ident, p.adjust.method = "bonferroni")


txnrd1 <- VlnPlot(combined, features = "TXNRD1", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  geom_boxplot(width=0.3) +
  theme(legend.position="none")
ggsave("new_figures/96h_TopBottom_TXNRD1_Vln.svg")


taldo1 <- VlnPlot(combined, features = "TALDO1", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  geom_boxplot(width=0.3) +
  theme(legend.position="none")
ggsave("new_figures/TopBottom_TALDO1_Vln.svg")

taldo1$data %>% dunn_test(TALDO1 ~ ident, p.adjust.method = "bonferroni")

sod1 <- VlnPlot(combined, features = "SOD1", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  geom_boxplot(width=0.3) +
  theme(legend.position="none")
ggsave("new_figures/TopBottom_SOD1_Vln.svg")

sod1$data %>% dunn_test(SOD1 ~ ident, p.adjust.method = "bonferroni")


PFKFB3 <- VlnPlot(combined, features = "PFKFB3", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  geom_boxplot(width=0.3) +
  theme(legend.position="none")
ggsave("new_figures/TopBottom_PFKFB3_Vln.svg")

PFKFB3$data %>% dunn_test(PFKFB3 ~ ident, p.adjust.method = "bonferroni")



library(escape)

GS.hallmark <- getGeneSets(library="H", species = "Homo sapiens")

combined <- runEscape(combined, 
                      method = "ssGSEA",
                      gene.sets = GS.hallmark, 
                      groups = 1000, 
                      min.size = 5,
                      new.assay.name = "escape.H")


ros <- VlnPlot(combined, features = "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY", assay="escape.H", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score') 
ggsave(ros, filename='new_figures/TopBottom_ROS_Vln.svg')

ros$data %>% dunn_test(`HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY` ~ ident, p.adjust.method = "bonferroni")


oxphos <- VlnPlot(combined, features = "HALLMARK-OXIDATIVE-PHOSPHORYLATION", assay="escape.H", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score') 
ggsave(oxphos, filename='new_figures/TopBottom_OXPHOS_Vln.svg')

oxphos$data %>% dunn_test(`HALLMARK-OXIDATIVE-PHOSPHORYLATION` ~ ident, p.adjust.method = "bonferroni")



glycolysis <- VlnPlot(combined, features = "HALLMARK-GLYCOLYSIS", assay="escape.H", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) + 
  geom_boxplot(width=0.3) + 
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score') 
ggsave(glycolysis, filename = "TopBottom_Glycolysis_Vln.svg")

glycolysis$data %>% dunn_test(`HALLMARK-GLYCOLYSIS` ~ ident, p.adjust.method = "bonferroni")

GS.Reactome <- getGeneSets(library="C2", subcategory = "CP:REACTOME", species = "Homo sapiens") 

combined <- runEscape(combined, 
                      method = "ssGSEA",
                      gene.sets = GS.Reactome, 
                      groups = 1000, 
                      min.size = 5,
                      new.assay.name = "escape.Reactome")


reac_ppp <- VlnPlot(combined, features = "REACTOME-PENTOSE-PHOSPHATE-PATHWAY", assay="escape.Reactome", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) + 
  geom_boxplot(width=0.3) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score')
reac_ppp$data %>% dunn_test(`REACTOME-PENTOSE-PHOSPHATE-PATHWAY` ~ ident, p.adjust.method = "bonferroni")


GS.Biocarta <- getGeneSets(library="C2", subcategory = "CP:BIOCARTA", species = "Homo sapiens") 
combined <- runEscape(combined, 
                      method = "ssGSEA",
                      gene.sets = GS.Biocarta, 
                      groups = 1000, 
                      min.size = 5,
                      new.assay.name = "escape.Biocarta")

arenrf2 <- VlnPlot(combined, features = "BIOCARTA-ARENRF2-PATHWAY", assay="escape.Biocarta", pt.size=0, cols = c("#ED1E2D", "#00AEEF", "#FF9999", "#99CCFF")) + 
  geom_boxplot(width=0.3) +
  stat_compare_means(method = "kruskal.test", label.x=2) +
  theme(legend.position="none") +
  xlab("") + ylab('Enrichment Score')
ggsave(arenrf2, filename = "TopBottom_ARENRF2_Vln.svg")

arenrf2$data %>% dunn_test(`BIOCARTA-ARENRF2-PATHWAY` ~ ident, p.adjust.method = "bonferroni")