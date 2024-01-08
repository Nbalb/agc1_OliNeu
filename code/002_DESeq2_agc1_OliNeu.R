# Following 
# https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://github.com/hbctraining/DGE_workshop_salmon_online
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# 
# Import salmon quant.sf files using tximeta
library(ComplexHeatmap)
library(EnhancedVolcano) 
library(ggbeeswarm) 
library(org.Mm.eg.db)
library(scales)
library(RColorBrewer)
library(DESeq2)
library(tximeta)
library(AnnotationHub)
library(apeglm)
library(tidyverse)
source("funcs/p_star.R")
dir.create("plots/")
dir <- "data/salmon"

files <- file.path(list.files(dir, "_quant", full.names = T), "quant.sf")
file.exists(files)

newNames <- data.frame("Sample" = str_remove(list.files(dir, "_quant"), "_quant"),
                       "genotype" = str_remove(list.files(dir, "_quant"), "[1-3]_quant"))

# Sample genotype
# 1    kd1       kd
# 2    kd2       kd
# 3    kd3       kd
# 4    wt1       wt
# 5    wt2       wt
# 6    wt3       wt

coldata <- data.frame(
  "names" = newNames$Sample,
  "files" = files)

fname <- "data/002_gse.rds"
if(!file.exists(fname)){
  dir.create("data/BFC")
  dir.create("data/AHC")
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
gse$gt <- factor(str_remove(gse$names, "[1-3]$")) |> 
  relevel("wt")
dds <- DESeqDataSet(gse, design = ~gt)
nrow(dds)   
# 35682

fname <- ("data/002_dds_rlog_normalized_agc1v2.rds")
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
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("plots/002_1_Euclidean_distance.png", h = 2000, w = 3000, res = 600)
Heatmap(sampleDistMatrix,
        col = colors, 
        row_labels = str_replace(rownames(sampleDistMatrix), "kd", "siAgc1 ") |> 
          str_replace("wt", "control "),
        column_title = "Overall similarity between samples",
        name = "Euclidean \ndistance\n",
        rect_gp = gpar(col = "grey60", lwd = 1)
)
dev.off()

sample_colors <- c(viridis_pal(option = "C", end = 0.8)(8)[1:3],
                   viridis_pal(option = "C", end = 0.8)(8)[6:8])
sample_labels <- c(paste0("siAgc1 rep", 1:3), paste0("control rep", 1:3))

plotPCA(rld, intgroup = c("names", "gt"), ntop = nrow(rld)) +
  labs(title = "Principal Component Analysis") +
  geom_point(size = 5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = sample_colors, 
                     labels = sample_labels,
                     name = "Sample") +
  theme_light()
ggsave("plots/002_2_PCA_plot.png", h = 1000, w = 1500, units = "px")

# Differential expression analysis
fname <- "data/002_des.rds"
if(!file.exists(fname)){
  des <- DESeq(dds)
  saveRDS(des, fname)
}else{des <- readRDS(fname)}

fname <- "data/002_res.rds"
if(!file.exists(fname)){
  
  res <- lfcShrink(des, coef = "gt_kd_vs_wt", type = "apeglm")
  
  if(all.equal(rowRanges(dds)$gene_id, rownames(res))){    # TRUE means that genes are all in the same order
    res$symbol <- rowRanges(dds)$symbol          # Add gene symbols
  }else{warning("genes in dds object and res object are not in the same order")}
  
  saveRDS(res, fname)
  res_tbl <- res |> 
    as.data.frame() |> 
    rownames_to_column("ensembl") |> 
    as_tibble()
  saveRDS(res_tbl,
          "data/002_res_tbl.rds")
  write_delim(res_tbl, "data/002_res_tbl.tsv", "\t")
  
}else{res <- readRDS(fname)}

# MA plot no shrinkage
png(paste0("plots/002_4_agc1_v2_MAplot_no_shrink.png"), w = 1600, h = 1200, res = 300)
plotMA(results(des), main = paste("MA plot siAgc1 OliNeu cells - no shrinkage"),
       alpha = 0.05,
       colSig = "#1E20FF80",
       colNonSig = "#88888870"
)
dev.off()

# MA plot shrinkage
png(paste0("plots/002_4_agc1_v2_MAplot.png"), w = 1600, h = 1200, res = 300)
plotMA(res, main = paste("MA plot siAgc1 OliNeu cells - apeglm shrinkage"),
       alpha = 0.05,
       colSig = "#1E20FF80",
       colNonSig = "#88888870",
       colLine = "#ff000080"
)
dev.off()

# Count number of significant genes
res |> as.data.frame() |> count(padj < 0.05)
# padj < 0.05     n
# 1       FALSE 12375
# 2        TRUE  2953
# 3          NA 20354

# Dispersion estimates
png("plots/002_5_Dispersion_estimates.png", w = 4000, h = 3000, res = 600)
plotDispEsts(des, main = "Dispersion estimates for siAgc1 OliNeu cells")
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
                             intgroup = c("gt", "names"), 
                             returnData = T) |>
      mutate(gene = factor(x))
    
  })  |>
  left_join(ord |> select(Ensembl_id, symbol), 
            by = c("gene" = "Ensembl_id")) |> 
  select(-gene) |> 
  rename(gene = symbol) |> 
  mutate(gene = factor(gene, levels = ord |> pull(symbol)),
         count = count + 0.5)


png(paste0("plots/002_6_agc1_top_genes_counts.png"), h = 3000, w = 6000, res = 600)

p <-  ggplot(best, aes(x = gt, y = count, color = names)) +
  scale_y_log10() +
  geom_beeswarm(size = 5, cex = 10, alpha = 0.8) +
  facet_wrap(~gene, nrow = 2) +
  labs(title = paste0("Best DE genes' counts in siAgc1 OliNeu cells"),
       x = "Genotype",
       y = parse(text=paste("log[10]","~normalized ~counts"))) +
  scale_x_discrete(labels = function(x) str_replace(x, "_", "\n")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = sample_colors)
print(p)

dev.off()

# Plot best scoring genes
ends <- res |>
  as.data.frame() |>
  rownames_to_column("ensembl_id") |>
  filter(padj < 0.05) |>
  filter(str_detect(ensembl_id, "^CG|^CR", negate = T)) |>    # Keep only known genes
  mutate(Sign = if_else(log2FoldChange > 0, "Up","Down") |>
           as.factor()) |> 
  group_by(Sign) |> 
  slice_max(abs(log2FoldChange), n = 25) |> 
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

png(paste0("plots/002_7_best_log2FC_genes.png"), h = 4000, w = 5500, res = 600)

p <- ggplot(ends, aes(x = c(1:50), y = log2FoldChange, fill = Sign)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Most significant differentially expressed genes in siAgc1 OliNeu cells", 
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
  coord_flip() +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  geom_text(aes(label = paste0(symbol, " - ", genedesc) |> str_wrap(55),
                y = ifelse(sign(log2FoldChange) > 0, -0.5, 0.5),
                hjust = ifelse(sign(log2FoldChange) > 0, 1, 0)),
            position = position_dodge(width = 0.2),
            size = 3.2,
            lineheight = 0.85)
print(p)

dev.off()


# Volcano plots
x <- res |> 
  as.data.frame() |> 
  rownames_to_column()

EnhancedVolcano(x, subtitle = "padj cutoff = 0.05 \nlog2FC cutoff = 1",
                lab = x$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0('Differential expression in siAgc1 OliNeu cells'),
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
                                filter(abs(log2FoldChange) > 5 | padj < 1e-150) |> 
                                pull(symbol))
) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))
ggsave("plots/002_8_Volcano_plot_agc1.png", h = 3500, w = 4500, units = "px", dpi = 600)

# Plot genes of interest
targets <-  c("Slc25a12", "Slc25a13", "Nat8l", "Aspa", "Acss1", "Srebf1",
              "Srebf2", "Fasn", "Olig1", "Olig2", "Olig3", "Sox10",
              "Pdgfra", "Cspg4", "Cldn11", "Mog", "Hdac9", "Hdac5", "Hdac2",
              "Hdac7", "Hdac3", "Hdac4", "Hdac1", "Hdac6", "Hdac11", "Hdac10", 
              "Hdac8", "Rai1")

fname <- "data/002_tpms_tbl.rds"
if(!file.exists(fname)){
  goi_full <- assays(des)[["abundance"]] |> 
    as.data.frame() |> 
    rownames_to_column("gene_id") |> 
    mutate(symbol = as.factor(rowRanges(des)$symbol),
           padj = res$padj,
           log2FC = res$log2FoldChange) |> 
    pivot_longer(starts_with("wt") | starts_with("kd"), names_to = "sample", values_to = "tpms")
  saveRDS(goi_full, fname)
}else{goi_full <- readRDS(fname)}
goi <- goi_full |> 
  filter(symbol %in% targets) |> 
  mutate(gt = str_remove_all(sample, "[0-9]"))

dir.create("plots/002_9_Genes_of_interest", showWarnings = F)

for(i in targets){
  goplot <- goi |> 
    filter(symbol == i) |> 
    mutate(gt = factor(gt, levels = c("wt", "kd"), labels = c("control", "siAgc1")))
  
  png(paste0("plots/002_9_Genes_of_interest/", i, ".png"), h = 2000, w = 2500, res = 600)
  p <- ggplot(goplot, aes(x = gt, 
                          y = tpms, 
                          color = sample)) +
    geom_beeswarm(size = 5, cex = 3) +
    labs(title = paste0(i, " expression \nin siAgc1 OliNeu cells"),
         subtitle = paste0(
           "log2FC = ",
           goplot |> pull(log2FC) |> head(1) |> scales::scientific(),
           "    padj = ",
           goplot |> pull(padj) |> head(1) |> scales::scientific()
         ),
         x = "Genotype", 
         y = "TPMs") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(labels = function(x) str_replace(x, "_", "\n")) +
    scale_color_manual(values = sample_colors,
                       labels = sample_labels,
                       name = "Sample") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = 
                                            ifelse(
                                              goplot |> 
                                                slice_head() |> 
                                                pull(padj) < 0.05,
                                              ifelse(goplot |> 
                                                       slice_head() |> 
                                                       pull(log2FC) < 0, 
                                                     "#8EDEF640", 
                                                     "#F8252915"),
                                              "grey90")))
  print(p)
  dev.off()
}

# Find mean TPMS for agc1 and agc2
goi_full |> 
  filter(str_detect(symbol, "Slc25a12|Slc25a13")) |> 
  mutate(pheno = str_remove(sample, "[0-9]")) |> 
  group_by(symbol, pheno) |> 
  summarise(mean = mean(tpms)) |> 
  ungroup()

# Plot heatmap of genes of interest
mtx <- goi |> 
  select(symbol, tpms, sample) |> 
  pivot_wider(names_from = symbol, values_from = tpms) |> 
  column_to_rownames("sample") |> 
  as.matrix() |> 
  scale() |> 
  t()
Heatmap(mtx, name = "TPMs", color_space = "rgb", column_names_rot = 0, 
        cluster_columns = F)

# Plot heatmap of best scoring genes
best_de <- res |> 
  as.data.frame() |> 
  na.omit() |> 
  mutate(sign = ifelse(log2FoldChange < 0, "neg", "pos")) |> 
  group_by(sign) |> 
  slice_max(abs(log2FoldChange), n = 25,) |> 
  ungroup() |> 
  pull(symbol)

mtx <- goi_full |> 
  select(symbol, tpms, sample) |> 
  filter(symbol %in% best_de) |> 
  pivot_wider(names_from = symbol, values_from = tpms) |> 
  column_to_rownames("sample") |> 
  as.matrix() |> 
  scale() |> 
  t()

logfc_bar <- res |> 
  as_tibble() |> 
  filter(symbol %in% best_de) |> 
  arrange(desc(log2FoldChange))

mtx <- mtx[logfc_bar$symbol,]

ha <- rowAnnotation(log2FC = anno_numeric(logfc_bar$log2FoldChange |> 
                                            round(2)), 
                    annotation_name_rot = 0)

png("plots/002_Heatmap_best_de.png", h = 2500, w = 3500, res = 350)
hm <- Heatmap(mtx, 
              cluster_columns = F, 
              cluster_rows = F,
              column_labels = str_replace(colnames(mtx), "wt", "control ") |> 
                str_replace("kd", "siAgc1 "),
              column_names_centered = T,
              column_names_rot = 0,
              column_title = "Top 30 differentially expressed genes",
              col = viridis_pal()(11),
              name = "scaled \nTPMs",
              right_annotation = ha,
              row_labels = paste0("  ", rownames(mtx)),
              row_names_centered = T
)
draw(hm)
dev.off()

# Phenotype heatmap
neuroep <- c("Nes", "Sox2", "Notch1", "Hes1", "Hes3", "Cdh1", "Ocln")
radial_glia <- c("Vim", "Nes", "Pax6", "Hes1", "Hes5", "Gfap", "Glast", "Blbp", "Tnc", "Cdh2", "Sox2") 
intermediate_progenitors <- c("Eomes", "Mash1", "Ascl1")
immature_neurons <- c("Dcx", "Tubb3", "Neurod1", "Tbr1", "Stmn1") # signifcant but unclear
opc <- c("Pdgfra", "Cspg4")
mature_oligodendrocytes <- c("Olig1", "Olig2", "Olig3", "Mbp", "Cldn11", "Mog", "Sox10")
schwann_cells <- c("Mpz", "Ncam1", "Gap43", "S100b", "Ngfr") # down
astrocytes <- c("Gfap", "Slc1a3", "Slc1a2", "Glul", "S100b", "Aldh1l1") #up
microglia <- c("Tmem119", "Itgam", "Ptprc", "Iba1", "Cx3cr1", "Adgre1", "Cd68", "Cd40")
mature_neurons <- c("Rbfox3", "Map2", "Nefm", "Nefh", "Syp", "Sypl1", "Sypl2", "Dlg4",
                    "Eno2", "Tubb3", "Nefl", "Gap43") #up
gaba_neurons <- c("Slc2a1", "Gabbr1", "Gabbr2", "Gad2", "Gad1") # 1 signif
glut_neurons <- c("Slc17a7", "Slc17a6", "Nmdar1", "Grin2b", "Gls", "Gls2", "Glul")

list_sig <- list(opc, mature_oligodendrocytes, schwann_cells, astrocytes, mature_neurons) |> 
  set_names("opc", "mature_oligodendrocytes", "schwann_cells", "astrocytes", "mature_neurons")
sig_tbl <- map_dfr(names(list_sig), function(x){
  res |> 
    as_tibble() |> 
    filter(symbol %in% !!sym(x), 
           !is.na(log2FoldChange),
           !is.na(padj)) |> 
    mutate(gs_name = x)
}) 

mtx <- goi_full |> 
  select(symbol, tpms, sample) |> 
  filter(symbol %in% sig_tbl$symbol) |> 
  pivot_wider(names_from = symbol, values_from = tpms) |> 
  column_to_rownames("sample") |> 
  as.matrix() |> 
  scale() |> 
  t()

mtx <- mtx[sig_tbl$symbol,]

ha <- rowAnnotation(space = anno_empty(border = F),
                    log2FC = anno_numeric(sig_tbl$log2FoldChange |> 
                                            round(2)),
                    Padj = anno_text(p_star(sig_tbl$padj), show_name = T),
                    annotation_name_rot = 0
                    )
rsplit <- sig_tbl |> 
  mutate(gs_name = str_replace_all(gs_name, "_", " ") |> 
           str_to_title() |> 
           str_wrap(10)) |> 
  pull(gs_name)

png("plots/002_Heatmap_best_brain_cell_markers.png", h = 2200, w = 3000, res = 350)
hm <- Heatmap(mtx, 
              cluster_columns = F, 
              cluster_rows = F,
              column_labels = str_replace(colnames(mtx), "wt", "control ") |> 
                str_replace("kd", "siAgc1 "),
              column_names_centered = T,
              column_names_rot = 0,
              column_title = "Brain cells markers expression",
              col = viridis_pal()(11),
              name = "scaled \nTPMs",
              right_annotation = ha,
              row_labels = paste0("  ", rownames(mtx)),
              row_names_centered = T,
              row_split = rsplit,
              row_gap = unit(2, "mm")
              
)
draw(hm)
dev.off()
