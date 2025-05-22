# scRNAseq data from:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6713459/
# Green CD, Ma Q, Manske GL, Shami AN, Zheng X, Marini S, Moritz L, Sultan C, Gurczynski SJ, Moore BB, Tallquist MD, Li JZ, Hammoud SS.
# A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell RNA-Seq. Dev Cell. 2018 Sep 10;46(5):651-667.e10. doi: 10.1016/j.devcel.2018.07.025. Epub 2018 Aug 23. PMID: 30146481; PMCID: PMC6713459.

# Load libraries
pacman::p_load(Seurat, tidyverse, pheatmap, cowplot, ggpubr, viridis, plyr, GenomicFeatures, txdbmaker, TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load count matrix
dir.create(path = "./Plots_PMC6713459", showWarnings = F)
count_matrix <- read.table(file = "./data_PMC6713459/GSE112393_MergedAdultMouseST25_DGE.txt", header = T)
data_PMC6713459.obj <- Seurat::CreateSeuratObject(counts = count_matrix, project = "Testis", min.cells = 3)
remove(count_matrix)

# Load metadata and add it to seurat object
attributes <- read_tsv(file = "./data_PMC6713459/GSE112393_MergedAdultMouseST25_PerCellAttributes.txt", col_names = T, skip = 3)
data_PMC6713459.obj@meta.data <- cbind(data_PMC6713459.obj@meta.data, attributes)

# Cell types from PMC6713459
cell_type <- c(
  "1" = "InnateLymph", "2" = "Macrophage", "3" = "Endothelial", "4" = "Myoid",
  "5" = "Leydig", "6" = "Sertoli", "7" = "Unknown",
  "8" = "SPG", "9" = "Scytes", "10" = "STids", "11" = "Elongating"
)

data_PMC6713459.obj@meta.data <- data_PMC6713459.obj@meta.data %>%
  mutate(cell_name = cell_type[data_PMC6713459.obj@meta.data$CellType] %>% as.vector())

VlnPlot(data_PMC6713459.obj, features = c("nFeature_RNA", "nCount_RNA", "%MT"), ncol = 3)

# 3. Normalize data ----------
data_PMC6713459.obj <- NormalizeData(data_PMC6713459.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Identify highly variable features --------------
data_PMC6713459.obj <- FindVariableFeatures(data_PMC6713459.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 100 most highly variable genes
top100 <- head(VariableFeatures(data_PMC6713459.obj), 100)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_PMC6713459.obj)
LabelPoints(plot = plot1, points = top100, repel = TRUE)

# 5. Scaling data-------------
all.genes <- rownames(data_PMC6713459.obj)
data_PMC6713459.obj <- ScaleData(data_PMC6713459.obj, features = all.genes)

# 6. Perform Linear dimensionality reduction --------------
data_PMC6713459.obj <- RunPCA(data_PMC6713459.obj, features = VariableFeatures(object = data_PMC6713459.obj))

# Visualize PCA results
DimHeatmap(data_PMC6713459.obj, dims = 1:20, cells = 500, balanced = TRUE)

# Determine dimensionality of the data
ElbowPlot(data_PMC6713459.obj)

# 7. Clustering ----------------------------
data_PMC6713459.obj <- FindNeighbors(data_PMC6713459.obj, dims = 1:20)

# Try different resolutions for clustering
data_PMC6713459.obj <- FindClusters(data_PMC6713459.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1, 2))
View(data_PMC6713459.obj@meta.data)

# Run UMAP
data_PMC6713459.obj <- RunUMAP(data_PMC6713459.obj, dims = 1:20)

# Label individual clusters using PMC6713459 cell classification ----------------------------
Idents(data_PMC6713459.obj) <- "cell_name"
p.clusters.res0.3 <- DimPlot(data_PMC6713459.obj, reduction = "umap", label = TRUE, repel = T)

# Prepare data for gene expression histogramsof Fmr1 and Htt for SPG, Scytes and STids ----------------------------

# 1. Subset seurat objects by cell types
spg.obj <- subset(x = data_PMC6713459.obj, cell_name == "SPG")
scytes.obj <- subset(data_PMC6713459.obj, cell_name == "Scytes")
sctids.obj <- subset(data_PMC6713459.obj, cell_name == "STids")
elon.obj <- subset(data_PMC6713459.obj, cell_name == "Elongating")

# 2. Extract normalized expression data
spg.exp <- as.data.frame(as.matrix(spg.obj@assays$RNA$data))
scytes.exp <- as.data.frame(as.matrix(scytes.obj@assays$RNA$data))
sctids.exp <- as.data.frame(as.matrix(sctids.obj@assays$RNA$data))
elon.exp <- as.data.frame(as.matrix(elon.obj@assays$RNA$data))

spg.exp.t <- as.data.frame(t(spg.exp[c("Fmr1", "Htt"), ]))
spg.exp <- pivot_longer(spg.exp.t, cols = c("Fmr1", "Htt"), names_to = "gene", values_to = "expression")

scytes.exp.t <- as.data.frame(t(scytes.exp[c("Fmr1", "Htt"), ]))
scytes.exp <- pivot_longer(scytes.exp.t, cols = c("Fmr1", "Htt"), names_to = "gene", values_to = "expression")

sctids.exp.t <- as.data.frame(t(sctids.exp[c("Fmr1", "Htt"), ]))
sctids.exp <- pivot_longer(sctids.exp.t, cols = c("Fmr1", "Htt"), names_to = "gene", values_to = "expression")

elon.exp.t <- as.data.frame(t(elon.exp[c("Fmr1", "Htt"), ]))
elon.exp <- pivot_longer(elon.exp.t, cols = c("Fmr1", "Htt"), names_to = "gene", values_to = "expression")


spg.exp$gene <- as.factor(spg.exp$gene)
scytes.exp$gene <- as.factor(scytes.exp$gene)
sctids.exp$gene <- as.factor(sctids.exp$gene)
elon.exp$gene <- as.factor(elon.exp$gene)

# Histograms with mean lines, including cells with 0 read counts for Htt or Fmr1 ----------------------------

# 1. SPGs
mu.spg.exp <- plyr::ddply(.data = subset(spg.exp, expression > -1), "gene", summarise, grp.mean = mean(expression))

p3 <- ggplot(data = subset(spg.exp, expression > -1)) +
  geom_histogram(aes(
    x = expression, y = ..count.. / sum(..count..),
    group = gene, fill = gene
  ), bins = 40, alpha = 0.5, position = "dodge") +
  geom_vline(data = mu.spg.exp, aes(xintercept = grp.mean), linetype = "dashed", color = c("#7AD151FF", "#414487FF")) +
  ylab(label = "Proportion of cells") +
  xlab(label = "Gene expression") +
  labs(title = "SPG (all cells)") +
  xlim(0.0000000001, max(spg.exp$expression)) + # ylim(0,0.017) +
  theme_classic() +
  scale_fill_manual(values = c("#7AD151FF", "#414487FF"), name = NULL)
print(p3)

ggsave(filename = "./Plots_PMC6713459/data_PMC6713459_spg_Fmr1_Htt_histogram_with_0.png", plot = p3, width = 4, height = 4)

# 2. Spermatocytes
mu.scytes.exp <- plyr::ddply(.data = subset(scytes.exp, expression > -1), "gene", summarise, grp.mean = mean(expression))

p4 <- ggplot(data = subset(scytes.exp, expression > -1)) +
  geom_histogram(aes(
    x = expression, y = ..count.. / sum(..count..),
    group = gene, fill = gene
  ), bins = 40, alpha = 0.5, position = "dodge") +
  geom_vline(data = mu.scytes.exp, aes(xintercept = grp.mean), linetype = "dashed", color = c("#7AD151FF", "#414487FF")) +
  ylab(label = "Proportion of cells") +
  xlab(label = "Gene expression") + # ylim(0,0.013) +
  labs(title = "Spermatocytes (all cells)") +
  xlim(0.0000000001, max(scytes.exp$expression)) +
  theme_classic() +
  scale_fill_manual(values = c("#7AD151FF", "#414487FF"), name = NULL)
print(p4)

ggsave(filename = "./Plots_PMC6713459/data_PMC6713459_Spermatocytes_Fmr1_Htt_histogram_with_0.png", plot = p4, width = 4, height = 4)


mu.sctids.exp <- plyr::ddply(.data = subset(sctids.exp, expression > -1), "gene", summarise, grp.mean = mean(expression))

# 3. Spermatids
p5 <- ggplot(data = subset(sctids.exp, expression > -1)) +
  geom_histogram(aes(
    x = expression, y = ..count.. / sum(..count..),
    group = gene, fill = gene
  ), bins = 40, alpha = 0.5, position = "dodge") +
  geom_vline(data = mu.sctids.exp, aes(xintercept = grp.mean), linetype = "dashed", color = c("#7AD151FF", "#414487FF")) +
  ylab(label = "Proportion of cells") +
  xlab(label = "Gene expression") +
  labs(title = "Spermatids (all cells)") +
  xlim(0.0000000001, max(sctids.exp$expression)) + # ylim(0,0.005) +
  theme_classic() +
  scale_fill_manual(values = c("#7AD151FF", "#414487FF"), name = NULL)
print(p5)

ggsave(filename = "./Plots_PMC6713459/data_PMC6713459_Spermatids_Fmr1_Htt_histogram_with_0.png", plot = p5, width = 4, height = 4)