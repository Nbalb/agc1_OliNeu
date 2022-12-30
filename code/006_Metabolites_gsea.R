library(msigdbr)
library(fgsea)
library(tidyverse)

# Import metabolite table
metab_orig <- read_delim("data/metabolites_quantitative_values_siagc1_vs_control.csv", delim = ";")
metab <- metab_orig |> 
  rename_with(~gsub("Err_", "Err.", .x), .cols = starts_with("Err")) |> 
  rename_with(~paste0("Concentration.", .x), .cols = starts_with(c("Cells", "Medium"))) |> 
  pivot_longer(cols = starts_with(c("Concentration", "Err")), 
               names_to = c(".value", "Fraction"), 
               names_sep = "\\.") |> 
  group_by(Metabolite) |> 
  mutate(group = case_when(
    mean(Concentration) < 20 ~ "0-20 ng",
    mean(Concentration) < 70 ~ "20-70 ng",
    mean(Concentration) > 70 ~ "70+ ng"
  )) |> 
  ungroup()

metab |> 
  filter(!(Metabolite %in% c("L-Glutamine", "ATP"))) |> 
  ggplot(aes(Metabolite, Concentration, fill = Fraction)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = Concentration-Err, ymax = Concentration+Err),
                position = position_dodge(.9),
                width = .5) +
  facet_wrap(~group, scales = "free") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#fcbcb8", "#97E5E7")) +
  labs(title = "Metabolites Levels (ng)")
ggsave("plots/010_Metabolites_concentration.png", h = 1500, w = 4000, units = "px")

# Check enrichment in Lasorsa gene sets
res <- readRDS("data/002_res.rds") |> 
  as_tibble()
stat <- res |> 
  na.omit() |> 
  filter(!duplicated(symbol)) |> 
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj),
         stat = -log10(padj) * log2FoldChange) |> 
  pull(stat, symbol)
  
pathways <- read_csv("results/pathways_to_search_Lasorsa.csv", col_names = c("pathway", "gene_name", "human_gene_symbol")) # Need to convert to mouse symbols
mmdb <- msigdbr(species = "Mus musculus")
pathways_mm <- pathways |> 
  inner_join(mmdb |> 
               select(human_gene_symbol, gene_symbol) |> 
               distinct(human_gene_symbol, .keep_all = T), by = "human_gene_symbol")
mlist <- split(pathways_mm$gene_symbol, pathways_mm$pathway)
fgsea_res <- fgseaMultilevel(pathways = mlist,
                             stats = stat, 
                             eps = 0)

# Since we didn't find any pathway to be significant we will just check the fold 
# changes for the genes of each pathway and 