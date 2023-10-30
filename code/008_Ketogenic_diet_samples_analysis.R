library(DESeq2)
library(tximeta)
library(AnnotationHub)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggbeeswarm)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(tidyverse)

# Read in raw counts
dir <- "data/008/salmon/"
files <- file.path(list.files(dir, "_quant", full.names = T), "quant.sf")
file.exists(files)

newNames <- data.frame("Sample" = str_remove(list.files(dir, "_quant"), "_quant"),
                       "diet" = str_remove(list.files(dir, "_quant"), "_R[1-5]_quant"))
newNames

coldata <- data.frame(
  "names" = newNames$Sample,
  "files" = files)

fname <- "data/008/001_gse.rds"
if(!file.exists(fname)){
  setTximetaBFC("data/BFC")
  setAnnotationHubOption("CACHE", "data/AHC")
  se <- tximeta(coldata)
  gse <- summarizeToGene(se)
  saveRDS(gse, file = fname)
}else(gse <- readRDS(fname))

round(colSums(assay(gse))/1e6, 1)
# kd1  kd2  kd3  wt1  wt2  wt3 
# 36.5 27.0 34.8 32.8 34.6 24.7 

# Data preparation and QC
set.seed(3)
gse$diet <- factor(str_remove(gse$names, "_R[1-5]$")) |> 
  relevel("Standard")
dds <- DESeqDataSet(gse, design = ~diet)
nrow(dds)   
# 35682

fname <- ("data/008/002_dds_rlog_normalized_agc1v2.rds")
if(!file.exists(fname)){
  rld <- rlog(dds, blind = FALSE)   # n<30 so we use rlog instead of vst
  saveRDS(rld, fname)
  all.equal(rowRanges(rld)$gene_id, rownames(rld))    # TRUE
  rownames(rld) <- rowRanges(rld)$gene_name
  saveRDS(assay(rld), "data/002_dds_rlog_normalized_assayonly.rds")
}else{rld <- readRDS(fname)}

sampleDists <- dist(t(assay(rld)))    # Euclidean distance between samples
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewerpal(9, "Blues")))(255)

png("plots/008/008_1_Euclidean_distance.png", h = 3000, w = 4200, res = 600)
Heatmap(sampleDistMatrix,
        col = colors,
        column_title = "Overall similarity between samples",
        name = "Euclidean \ndistance\n",
        rect_gp = gpar(col = "grey60", lwd = 1)
)
dev.off()

sample_colors <- c(brewer.pal(7, "Blues")[3:7], brewer.pal(7, "Greens")[3:7])

plotPCA(rld, intgroup = c("names", "diet"), ntop = nrow(rld)) +
  labs(title = "Principal Component Analysis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = sample_colors, name = "Group", labels = rld$names) +
  theme_bw()
ggsave("plots/008/008_2_PCA_plot.png", h = 1500, w = 3000, units = "px")

# Differential expression analysis
fname <- "data/008/003_des.rds"
if(!file.exists(fname)){
  des <- DESeq(dds)
  saveRDS(des, fname)
}else{des <- readRDS(fname)}

fname <- "data/008/004_res.rds"
if(!file.exists(fname)){
  
  res <- lfcShrink(des, coef = "diet_Ketogenic_vs_Standard", type = "apeglm")
  
  if(all.equal(rowRanges(dds)$gene_id, rownames(res))){    # TRUE means that genes are all in the same order
    res$symbol <- rowRanges(dds)$symbol          # Add gene symbols
  }else{warning("genes in dds object and res object are not in the same order")}
  
  saveRDS(res, fname)
  res_tbl <- res |> 
    as.data.frame() |> 
    rownames_to_column("ensembl") |> 
    as_tibble()
  saveRDS(res_tbl, "data/008/002_res_tbl.rds")
  write_delim(res_tbl, "data/008/002_res_tbl.tsv", "\t")
  
}else{res <- readRDS(fname)}

# MA plot no shrinkage
png("plots/008/008_3_agc1_v2_MAplot_no_shrink.png", w = 2000, h = 1200, res = 300)
plotMA(results(des), main = paste("MA plot keto diet oligodendrocyte - no shrinkage"),
       alpha = 0.05,
       colSig = "#1E20FF80",
       colNonSig = "#88888870"
)
dev.off()

# MA plot shrinkage
png("plots/008/008_4_agc1_v2_MAplot.png", w = 2000, h = 1200, res = 300)
plotMA(res, main = paste("MA plot keto diet oligodendrocyte - apeglm shrinkage"),
       alpha = 0.05,
       colSig = "#1E20FF80",
       colNonSig = "#88888870",
       colLine = "#ff000080"
)
dev.off()

# Count number of significant genes
res |> as.data.frame() |> count(padj < 0.05)
# padj < 0.05     n
# 1       FALSE 12030
# 2        TRUE   149
# 3          NA 23503

# Dispersion estimates
png("plots/008/008_5_Dispersion_estimates.png", w = 4000, h = 3000, res = 600)
plotDispEsts(des, main = "Dispersion estimates for keto diet oligodendrocytes")
dev.off()

# Plot best genes' counts
ord <- res |> 
  as.data.frame() |> 
  rownames_to_column("Ensembl_id") |> 
  filter(padj < 0.05) |> 
  arrange(-log2FoldChange) |> 
  slice(-((6):(n()-5)))

best <- ord |>
  pull(Ensembl_id) |> 
  map_dfr(function(x){
    
    geneCounts <- plotCounts(dds = des, 
                             gene = x, 
                             intgroup = c("diet", "names"), 
                             returnData = T) |>
      mutate(gene = factor(x))
    
  })  |>
  left_join(ord |> select(Ensembl_id, symbol), 
            by = c("gene" = "Ensembl_id")) |> 
  select(-gene) |> 
  rename(gene = symbol) |> 
  mutate(gene = factor(gene, levels = ord |> pull(symbol)),
         count = count + 0.5)

ggplot(best, aes(x = diet, y = count, color = names)) +
  scale_y_log10() +
  ggbeeswarm::geom_beeswarm(size = 5, cex = 6, alpha = 0.8) +
  facet_wrap(~gene, nrow = 2) +
  labs(title = paste0("Best DE genes' counts in keto diet oligodendrocytes"),
       x = "Genotype",
       y = parse(text=paste("log[10]","~normalized ~counts"))) +
  scale_x_discrete(labels = function(x) str_replace(x, "_", "\n")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = sample_colors) +
  theme_bw()
ggsave("plots/008/008_6_agc1_top_genes_counts.png", h = 1500, w = 3000, units = "px")

# Plot best scoring genes
source("funcs/p_star.R")
n <- 10
ends <- res |>
  as.data.frame() |>
  rownames_to_column("ensembl_id") |>
  filter(padj < 0.05) |>
  mutate(Sign = if_else(log2FoldChange > 0, "Up","Down") |>
           as.factor()) |> 
  group_by(Sign) |> 
  slice_max(abs(log2FoldChange), n = n) |> 
  ungroup() |> 
  arrange(log2FoldChange)

gene_desc <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keys = ends$ensembl_id, 
                                   columns = c("ENSEMBL", "GENENAME"), 
                                   keytype = "ENSEMBL")

ends <- ends |> 
  left_join(gene_desc, by = c("ensembl_id" = "ENSEMBL")) |> 
  mutate(genedesc = str_replace(GENENAME, "(.*,.*,.*,.*),.*", "\\1"),
         GENENAME = NULL)

png(paste0("plots/008/002_7_best_log2FC_genes.png"), h = 4000, w = 5500, res = 600)

ggplot(ends, aes(x = seq(n*2), y = log2FoldChange, fill = Sign)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Most significant differentially expressed genes in keto diet oligodendrocytes", 
       subtitle = paste0("padj < ", signif(max(ends$padj), 4), " - only known genes are shown"),
       y = "log2FoldChange", 
       x = "") +
  geom_text(aes(label = p_star(padj)), 
            hjust = ifelse(ends$Sign == "Up", -0.3, 1.3),
            vjust = 0.75) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank()) +
  coord_flip(ylim = c(-3.5, 7.5)) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  geom_text(aes(label = paste0(symbol, " - ", genedesc) |> str_wrap(55),
                y = ifelse(sign(log2FoldChange) > 0, -0.1, 0.1),
                hjust = ifelse(sign(log2FoldChange) > 0, 1, 0)),
            position = position_dodge(width = 0.2),
            size = 3.2,
            lineheight = 0.85)
ggsave("plots/008/008_7_best_log2FC_genes.png", h = 2000, w = 3000, units = "px")

# Volcano plots
x <- res |> as.data.frame() |> rownames_to_column()

png(paste0("plots/008/008_8_Volcano_plot_ket_diet.png"), h = 3500, w = 4500, res = 600)
EnhancedVolcano(x, subtitle = "padj cutoff = 0.05 \nlog2FC cutoff = 1",
                lab = x$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0('Differential expression in keto diet oligodendrocytes'),
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 1,
                axisLabSize = 14,
                titleLabSize = 14,
                subtitleLabSize = 10,
                captionLabSize = 10,
                pointSize = 2.5,
                labSize = 3.5,
                legendLabSize = 10,
                drawConnectors = T,
                min.segment.length = 0,
                widthConnectors = 0.8,
                arrowheads = F,
                caption = paste0("Significantly downregulated genes = ", 
                                 x |> 
                                   filter(log2FoldChange < -1 & padj < 0.05) |>
                                   nrow(),
                                 "     ",
                                 "Significantly upregulated genes = ",
                                 x |> 
                                   filter(log2FoldChange > 1 & padj < 0.05) |> 
                                   nrow()),
                selectLab = c("Srebf1", "Srebf2", "Slc25a12", "Slc25a13",
                              x |> 
                                filter(log2FoldChange > 3 | padj < 1e-5 | log2FoldChange < -1.2) |> 
                                pull(symbol))
) +
  # coord_cartesian(xlim = c(-2.5, 7.5)) + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

dev.off()

