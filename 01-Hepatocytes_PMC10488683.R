# Bulk RNAseq data from hepatocytes from PMC10488683

# Load required libraries --------------
pacman::p_load(GenomicFeatures, Seurat, tidyverse, pheatmap, cowplot, ggpubr, viridis, plyr, org.Mmu.eg.db, edgeR)

# Calculate transcrip lengths from mmu10 (GRCm38.102) annotation for TPM normalization --------------
txdb <- txdbmaker::makeTxDbFromGFF("./data_PMC10488683/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff_NO_RETAIN_INTRON.gtf",
  format = "gtf"
)
gtf <- read_tsv(file = "./data_PMC10488683/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff_NO_RETAIN_INTRON.gtf", comment = "#", col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "description"))

gene.gtf <- subset(gtf, feature == "gene")
geneid2name <- gene.gtf %>%
  mutate(description = str_remove_all(string = description, pattern = "\"")) %>%
  separate(
    col = description, remove = T,
    sep = " |; ", into = c(as.character(1:19))
  ) %>%
  dplyr::select(c("2", "6")) %>%
  rename(all_of(c("2" = "gene_id", "6" = "gene_name"))) %>%
  mutate(gene_id = str_remove_all(string = gene_id, pattern = "\\.[0-9]+$"))

geneid2name <- read_tsv(file = "./data_PMC10488683/mart_export.txt", col_names = T)

exons.list.per.gene <- exonsBy(txdb, by = "gene")
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))
colnames(exonic.gene.sizes) <- c("length")
rownames(exonic.gene.sizes) <- rownames(exonic.gene.sizes) %>% str_remove_all(pattern = "\\.[0-9]+")
exonic.gene.sizes <- rownames_to_column(exonic.gene.sizes, var = "gene_id")
exonic.gene.sizes <- unique(merge(x = exonic.gene.sizes, y = geneid2name, by = "gene_id"))

# Download data from PMC10488683, file "fig4_7.csv" --------------
PMC10488683.raw <- read.delim(file = "https://digital.lib.washington.edu/researchworks/bitstreams/74b92b66-3445-4bbf-b2f6-db714bcf0241/download", header = T, sep = ",")

# Keep only WT samples and genes that we have transcript length info
PMC10488683.sub <- subset(x = PMC10488683.raw, gene %in% exonic.gene.sizes$gene_name)
PMC10488683.sub <- PMC10488683.sub %>% dplyr::select(starts_with(c("gene", "WT_")))

# Add transcript length
PMC10488683.sub <- merge(x = PMC10488683.sub, y = exonic.gene.sizes, by.x = "gene", by.y = "gene_name", all.x = T)

# Remove duplicated rows
PMC10488683.sub <- PMC10488683.sub[!duplicated(PMC10488683.sub), ]

# Remove solute carrier family 6  gene that is annotated as Htt with ID=ENSMUSG00000020838
PMC10488683.sub <- subset(PMC10488683.sub, !(gene == "Htt" & gene_id == "ENSMUSG00000020838"))

# Normalize to FPKM
PMC10488683.sub.fpkm <- edgeR::rpkm(
  y = PMC10488683.sub %>% dplyr::select(starts_with(c("WT_"))),
  gene.length = PMC10488683.sub$length, normalized.lib.sizes = F
) %>%
  as.data.frame()
PMC10488683.sub.fpkm$gene_name <- PMC10488683.sub$gene

# FPKM to TPM
PMC10488683.sub.tpm <- apply(dplyr::select(PMC10488683.sub.fpkm, starts_with("WT_")), 2, function(x) x * 10^6 / sum(x))
rownames(PMC10488683.sub.tpm) <- PMC10488683.sub$gene

# Extract expression for Htt and Fmr1
genes_of_interest.tpm <- subset(PMC10488683.sub.tpm, rownames(PMC10488683.sub.tpm) %in% c("Htt", "Fmr1"))

# Reshape genes_of_interest.tpm for plotting
genes_of_interest.tpm <- t(genes_of_interest.tpm) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(names_to = "gene", values_to = "tpm", cols = !sample)

genes_of_interest.tpm$gene <- factor(genes_of_interest.tpm$gene, levels = c("Fmr1", "Htt"))


# Running t.test paired statistics because N < 10 --------------

# Plot expressions
dir.create(path = "./Plots_PMC10488683", showWarnings = F)

p.PMC10488683 <- ggboxplot(genes_of_interest.tpm,
  x = "gene", y = "tpm",
  color = "gene", palette = c("#7AD151FF", "#414487FF"),
  add = "jitter", add.params = list(size = 1, alpha = 1),
  title = "Hepatocytes (PMC10488683)"
) +
  stat_compare_means(
    label.y = max(genes_of_interest.tpm$tpm) * 1.2,
    comparisons = list(c("Fmr1", "Htt")),
    paired = T, method = "t.test"
  ) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = max(genes_of_interest.tpm$tpm) * 1.1,
        label = paste("Mean=", round(mean(x), digits = 2))
      )
    },
    geom = "text"
  ) +
  theme(legend.position = "none") + ylab("TPM")

ggsave(filename = "./Plots_PMC10488683/PMC10488683_liver_Fmr1_Htt_boxplot.png", plot = p.PMC10488683, height = 4, width = 3)

print(p.PMC10488683)
