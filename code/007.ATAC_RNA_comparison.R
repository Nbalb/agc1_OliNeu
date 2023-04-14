library(ggrepel)
library(ggpointdensity)
library(tidyverse)

atac_de <- read_csv("../agc1v2/data/ATACSeq/005.1c_ATAC_de.csv")
dge <- readRDS("data/002_res_tbl.rds")

point_plot <- atac_de |> 
  filter(padj < 0.05) |> 
  inner_join(dge |> 
               filter(padj < 0.05), 
             by = c("gene_symbol" = "symbol"), suffix = c("_atac", "_rna")) |> 
  mutate(labels = ifelse(abs(log2FoldChange_atac) > 2 | abs(log2FoldChange_rna) > 2, 
                         gene_symbol, ""))
point_plot |> 
  ggplot(aes(log2FoldChange_rna, log2FoldChange_atac, label = labels)) +
  geom_pointdensity() +
  labs(title = "Significant genes comparison",
       x = "RNASeq log2FoldChange",
       y = "ATACSeq log2FoldChange",
       subtitle = paste0(nrow(point_plot), " genes compared")) + 
  scale_color_viridis_c(option = "C") +
  theme_bw()
ggsave("plots/007/001_logFC_ATAC_vs_RNA.png", h = 1200, w = 2000, units = "px")

# Density plot 
dens_plot <- point_plot |> 
  select(gene_symbol, starts_with("Log")) |> 
  pivot_longer(-gene_symbol, names_to = "experiment", names_prefix = "log2FoldChange_", values_to = "log2FC") |> 
  mutate(experiment = str_to_upper(experiment))

mw_res <- wilcox.test(dens_plot |> filter(experiment == "ATAC") |> pull(log2FC),
                      dens_plot |> filter(experiment == "RNA") |> pull(log2FC),
                      alternative = "two.sided")

dens_plot |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison.png", h = 1200, w = 2000, units = "px")

mw_res_pos <- wilcox.test(dens_plot |> filter(experiment == "ATAC", log2FC > 0) |> pull(log2FC),
                          dens_plot |> filter(experiment == "RNA", log2FC > 0) |> pull(log2FC),
                          alternative = "two.sided")

dens_plot |> 
  filter(log2FC > 0) |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of positive log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res_pos$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison_positive.png", h = 1200, w = 2000, units = "px")

mw_res_neg <- wilcox.test(dens_plot |> filter(experiment == "ATAC", log2FC < 0) |> pull(log2FC),
                          dens_plot |> filter(experiment == "RNA", log2FC < 0) |> pull(log2FC),
                          alternative = "two.sided")

dens_plot |> 
  filter(log2FC < 0) |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of negative log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res_neg$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison_negative.png", h = 1200, w = 2000, units = "px")
