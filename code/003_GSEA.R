library(fgsea)
library(data.table)
library(ComplexHeatmap)
library(scales)
library(DESeq2)
library(msigdbr)
library(tidyverse)
source("funcs/plotEnrichment2.R")
source("funcs/pretty_path_label.R")
source("funcs/p_star.R")

# Following msigdbr and fgsea vignettes ----
fname <- "data/003_mdf_mm.rds"
if(!file.exists(fname)){
  library(msigdbr) 
  msigdbr_collections() |> print(n = Inf)
  mdf <- msigdbr(species = "Mus musculus") |>     # Retrieve all Dm gene sets
    filter(gs_cat %in% c("C2", "C5", "H"), 
           gs_subcat %in% c("", "GO:MF", "GO:BP", "GO:CC", "CP:KEGG", "CP:WIKIPATHWAYS", "CP:REACTOME"))
  saveRDS(mdf, file = fname)
}else(mdf <- readRDS(fname))

mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

# Get signatures
res <- readRDS("data/002_res.rds")
fname <- "data/003_fgsea_results.rds"

if(!file.exists(fname)){
  
  df <- res |> as.data.frame() |> 
    filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(symbol)) |> 
    mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
  
  #signature <- setNames(-log10(df$padj)*sign(df$log2FoldChange), rownames(df))
  sig <- setNames(df$log2FoldChange, df$symbol)
  #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
  fgseaRes <- fgseaMultilevel(pathways = mlist,
                              stats = sig,
                              eps = 0,
                              nproc = parallel::detectCores()-2,
  )  
  
  
  saveRDS(fgseaRes, fname)
  
}else{fgseaRes <- readRDS(fname)}

# Plots ----
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      mlist, sig)

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

png(paste0("plots/003_01_agc1_collapsed_pathways.png"), h = 5000, w = 6500, res = 600)
plotGseaTable(mlist[mainPathways], sig, fgseaRes, 
              gseaParam = 0.5, colwidths = c(12,4,1.8,2,2))
dev.off()

## Barplot gsea
ends <- fgseaRes |>
  filter(pathway %in% mainPathways, NES > 0, padj < 0.05) |> 
  slice_min(padj, n = 10) |> 
  arrange(-NES) |> 
  bind_rows(fgseaRes |>
              filter(NES < 0 & padj < 0.05) |> 
              slice_min(padj, n = 10) |> 
              arrange(-NES)) |> 
  arrange(NES)

png(paste0("plots/003_02_agc1_GSEA_barplot.png"), h = 3000, w = 4000, res = 600)

plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES)))) +
  geom_bar(stat='identity') +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Best scoring pathways in siAgc1 OliNeu cells", y = "Combined Score", x="") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  coord_flip(ylim = c(-0.8, 2.1)) +
  geom_text(aes(label = pretty_path_label(pathway, wrap = 30), 
                y = ifelse(NES < 0, 0.05, -0.05),
                hjust = ifelse(NES < 0, 0, 1)),
            position = position_dodge(width = 0),
            size = 3,
            lineheight = 0.85) +
  geom_text(aes(label = p_star(padj)), 
            hjust = ifelse(sign(ends$NES) < 0, 1.3, -0.3),
            vjust = 0.75)


print(plt)
dev.off()

## Plot single GSEAs
pdf(paste0("plots/003_03_agc1_main_pathways_enrichment.pdf"), w = 5, h = 3)
for(o in 1:length(collPathways)){
  
  NES <- signif(fgseaRes[pathway == mainPathways[o]]$NES, 3)
  padj <- signif(fgseaRes[pathway == mainPathways[o]]$padj, 3)
  title <- paste0("GSEA in \n", 
                  pretty_path_label(mainPathways[o]) |> str_wrap(50))
  
  p <- plotEnrichment2(mlist[[mainPathways[o]]], stats = sig) +
    labs(title = title, 
         subtitle = paste0("NES = ", NES, "  p.adj = ", padj)) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5)) +
    xlim(0, 10500)
  print(p)
}
dev.off()

# Heatmap genes involved in myelination
# Check if pathways related to myelination are significant
myelin_pw <- fgseaRes |> 
  as_tibble() |> 
  filter(str_detect(pathway, "MYELIN")) |> 
  arrange(padj)

## Since none of the pathways are significant, let's check at gene level
myelin_genes <- mdf |> 
  filter(str_detect(gs_name, "MYELIN")) |> 
  pull(gene_symbol) |> 
  c("Srebf1", "Nat8l", "Aspa", "Acss", "Fasn", "Rai1")

res_myelin <- res |> 
  as_tibble(rownames = "ensembl_id") |> 
  filter(symbol %in% myelin_genes,
         padj < 0.05)

logfc_bar <- res |> 
  as_tibble() |> 
  filter(symbol %in% res_myelin$symbol) |> 
  mutate(p_star = p_star(padj)) |> 
  arrange(desc(log2FoldChange))

tpms <- readRDS("data/002_tpms_tbl.rds")
mtx <- tpms |> 
  select(symbol, tpms, sample) |> 
  filter(symbol %in% res_myelin$symbol) |> 
  pivot_wider(names_from = symbol, values_from = tpms) |> 
  column_to_rownames("sample") |> 
  as.matrix() |> 
  scale() |> 
  t()
mtx <- mtx[logfc_bar$symbol,]

ha <- rowAnnotation(log2FC = anno_numeric(logfc_bar$log2FoldChange |> 
                                             round(2)), 
                    Padj = anno_text(logfc_bar$p_star, show_name = T),
                    annotation_name_rot = 0)

png("plots/003_Heatmap_myelination_tpms.png", h = 2700, w = 3200, res = 350)
hm <- Heatmap(mtx, 
        cluster_columns = F, 
        cluster_rows = F,
        column_names_rot = 0,
        column_title = "Genes involved in myelination",
        col = viridis_pal()(11),
        name = "scaled \nTPMs",
        right_annotation = ha,
        row_labels = paste0("  ", rownames(mtx)),
        row_names_centered = T
        )
draw(hm)
dev.off()

### ClusterProfiler

### Msigdb
#### GSEA 
sig <- setNames(as.numeric(df_id$log2FoldChange), as.character(df_id$ENTREZID))
sig <- sort(sig, decreasing = T)

nmdf <- mdf |>
  dplyr::select(gs_name, entrez_gene) |> 
  dplyr::rename(term = gs_name, gene = entrez_gene) |> 
  as.data.frame()

gsemdf <- GSEA(sig, TERM2GENE = nmdf, eps = 0,)

gsemdf@result <- gsemdf@result |> 
  mutate(Description = pretty_path_label(Description, 50))

png("plots/003_04_Cluster_Profiler_msigdb_GSEA_dotplot.png", h = 4000, w = 6000, res = 600)
dotplot(gsemdf, showCategory = 15) + 
  labs(title = "MSigDB subset GSEA for siAgc1 OliNeu Cells") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Enricher
library(enrichR)

listEnrichrDbs()
dbs <- c("CCLE_Proteomics_2020", "ProteomicsDB_2020", "DisGeNET", "Pfam_Domains_2019", "KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "GWAS_Catalog_2019", "Rare_Diseases_AutoRIF_Gene_Lists", "Rare_Diseases_GeneRIF_Gene_Lists",
         "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "ARCHS4_TFs_Coexp")

upreg <- res |> 
  as.data.frame() |> 
  filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange), log2FoldChange > 0) |> 
  pull(symbol)

downreg <- res |> 
  as.data.frame() |> 
  filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange), log2FoldChange < 0) |> 
  pull(symbol)

eup <- enrichr(upreg, dbs)
edn <- enrichr(downreg, dbs)

dir.create("plots/003_Enrichr/")
for(x in dbs){
  
  pl <- eup[[x]] |> 
    filter(Adjusted.P.value < 0.05) |> 
    slice_max(Combined.Score, n = 10) |> 
    bind_rows(edn[[x]] |> 
                filter(Adjusted.P.value < 0.05) |> 
                slice_max(Combined.Score, n = 10) |> 
                mutate_at(vars(Combined.Score), ~ -(.)) |> 
                arrange(Combined.Score)) |> 
    mutate(Generatio = map_dbl(Overlap, 
                               ~eval(parse(text = .x), list(iter = 1))), 
           .after = Overlap) |> 
    mutate(Term = str_remove_all(Term, "human.*|TenP.*|\\(GO.*|WP.*|BTO.*") |> 
             str_wrap(40)) |> 
    dplyr::select(-Genes)
  
  if(nrow(pl) < 3) {next}
  
  ti <- str_replace_all(x, "_", " ")
  
  ggplot(pl, aes(reorder(Term, Combined.Score), Combined.Score, color = Generatio)) +
    geom_point(aes(size = -log10(Adjusted.P.value))) +
    geom_segment(aes(x = Term, xend = Term, y = 0, yend = Combined.Score)) +
    coord_flip() +
    scale_color_viridis_c(direction = -1) +
    theme_bw() +
    labs(title = paste0("Enrichment analysis for siAgc1 OliNeu Cells in \n", ti, " dataset"),
         size = expression(-log[10]~padj)) +
    xlab("Term") +
    ylab("Enrichr Combined Score") +
    scale_size_continuous(range = c(3,7.5)) +
    theme(axis.text.y = element_text(size = 12))
  ggsave(paste0("plots/003_Enrichr/003_Agc1_enrichment_", x, ".png"), width = 10, h = 6, dpi = 700)
  
}

