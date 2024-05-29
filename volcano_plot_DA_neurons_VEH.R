
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(mclust)


# Data loading and prep ---------------------------------------------------

se <- readRDS("data/SummExp_1191_CNS_CSS_VTA_neurons_VEH_Pryce.rds")

# Identify expressed and non-expressed genes via Gaussian mixture model on
# the median log2 counts  across samples
median_across_samples <- apply(log2(assay(se, "counts") + 1), 1, median)
g_mm <- Mclust(median_across_samples, G = 2, verbose = FALSE)
keep <- which(g_mm$classification == 2)

se_expr <- se[keep, ]

# Specify the reference level
se_expr$MFGroup <- factor(se_expr$MFGroup)
se_expr$MFGroup <- relevel(se_expr$MFGroup, "control_vehicle_DA_neuron")

# DESeq object construction
mtx <- as.matrix(assay(se_expr, "counts"))
dds <- DESeqDataSetFromMatrix(countData = round(mtx), 
                              colData = colData(se_expr), 
                              design = ~ MFGroup)

# DGEA
dds <- DESeq(dds)

# Results for specific contrasts
res <- results(dds, contrast=c("MFGroup", "CSS_vehicle_DA_neuron", "control_vehicle_DA_neuron"))
ttl <- "CSS Vehicle vs Control Vehicle"

# Check stats
sum(res$log2FoldChange > 0.5 & res$pvalue < 0.001, na.rm = TRUE)
sum(res$log2FoldChange < -0.5 & res$pvalue < 0.001, na.rm = TRUE)


# Volcano plots -----------------------------------------------------------

# Tibble for volcano plot (|log2FC > 0.5| and nominal p < 0.001)
for_volc <- as_tibble(res, rownames = "ensg") %>% 
  left_join(as_tibble(rowData(se_expr)), by = "ensg") %>% 
  dplyr::select(ensg, symbol, baseMean, log2FoldChange, pvalue, padj) %>%
  dplyr::rename(logFC = log2FoldChange,
                avgExpr = baseMean) %>% 
  mutate(status = case_when(
    logFC > 0.5 & pvalue < 0.001 ~ "Up-regulated",
    logFC < -0.5 & pvalue < 0.001 ~ "Down-regulated",
    TRUE ~ "Not significant"),
    lbl = case_when(
      abs(logFC) > 0.5 & pvalue < 0.001 ~ symbol,
      TRUE ~ ""))

# Volcano - labels
for_volc <- for_volc |> 
  mutate(lbl = ifelse(ensg == "ENSMUSG00000110038", "Gm45570", lbl))

# Volcano plot
vlbl <- ggplot(for_volc, aes(logFC, -log10(pvalue), colour = status, label = lbl)) +
  geom_point() +
  geom_text_repel(max.overlaps = 20) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "black") +
  geom_vline(xintercept = -0.5, linetype = 2, colour = "black") +
  geom_hline(yintercept = -log10(0.001), linetype = 2, colour = "black") +
  scale_color_manual(values = c("Up-regulated" = "red", 
                                "Down-regulated" = "blue", 
                                "Not significant" = "grey")) +
  labs(x = "log2 fold-change", y = "-log10(p-value)", 
       title = ttl) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0("figs/", ttl, "_VEH_7x6in.svg"), plot = vlbl, 
       width = 7, height = 6, units = "in")

