library(ComplexHeatmap)
library(msigdbr)
library(fgsea)
library(tidyverse)

# Load data
res_diet <- readRDS("data/008/002_res_tbl.rds")
res_kd <- readRDS("data/002_res_tbl.rds")

res <- res_kd |> 
  inner_join(res_diet, by = "ensembl", suffix = c("_kd", "_diet"))

res |> 
  filter(padj_kd < 0.05, padj_diet < 0.05) |> 
  mutate(diff = log2FoldChange_diet - log2FoldChange_kd) |> 
  select(matches("symbol|log2FoldChange|padj|diff")) |> 
  pivot_longer(!diff, names_to = c(".value", "contrast"), names_sep = "_") |> 
  ggplot(aes(reorder(symbol, -diff), log2FoldChange, fill = contrast)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "LFC differences in diet contrast vs kdAgc1 contrast",
       x = "Gene symbol") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/009/001_LFC_comparison.png", h = 1200, w = 2400, units = "px")

# Run GSEA on diet dataset
mdf <- readRDS("data/003_mdf_mm.rds")
mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

stat <- res_diet |> 
  na.omit() |> 
  filter(!duplicated(symbol)) |> 
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj),
         stat = -log10(padj) * log2FoldChange) |> 
  pull(stat, symbol)

fgseaRes <- fgseaMultilevel(pathways = mlist,
                            stats = stat,
                            eps = 0)  # No significant results

# Heatmap of most differentially expressed genes
gse_diet <- readRDS("data/008/001_gse.rds")
gse_kd <- readRDS("data/002_gse.rds")

gse_join <- bind_cols(rownames(gse_diet),
                      assays(gse_diet)$abundance, 
                      assays(gse_kd)$abundance)
