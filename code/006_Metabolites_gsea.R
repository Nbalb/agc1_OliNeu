library(msigdbr)
library(fgsea)
library(ggpattern)
library(tidyverse)

# Import metabolite table
metab_orig <- read_delim("data/metabolites_quantitative_values_siagc1_vs_control.csv", delim = ";")
metab <- metab_orig |> 
  rename_with(~gsub("Err_", "Err.", .x), .cols = starts_with("Err")) |> 
  rename_with(~paste0("Concentration.", .x), .cols = starts_with(c("Cells", "Medium"))) |> 
  pivot_longer(cols = starts_with(c("Concentration", "Err")), 
               names_to = c(".value", "Fraction"), 
               names_sep = "\\.") |> 
  mutate(Metabolite = ifelse(Metabolite == "N-Acetylaspartate", "NAA", Metabolite)) |> 
  group_by(Metabolite) |> 
  mutate(group = case_when(
    mean(Concentration) < 20 ~ "0-20 ng",
    mean(Concentration) < 70 ~ "20-70 ng",
    mean(Concentration) < 150 ~ "70-150 ng",
    mean(Concentration) > 150 ~ "150+ ng"
  ) |> 
    factor(levels = c("0-20 ng", "20-70 ng", "70-150 ng", "150+ ng"))) |> 
  ungroup() |> 
  separate_wider_delim(Fraction, names = c("Fraction", "Genotype"), delim = "_") |> 
  mutate(Genotype = ifelse(Genotype == "wt", "control", "siAgc1"))

metab |> 
  ggplot(aes(Metabolite, Concentration, fill = Genotype, pattern = Fraction)) +
  geom_bar_pattern(stat = "identity", 
                   position = position_dodge(), 
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 30,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin = Concentration-Err, ymax = Concentration+Err),
                position = position_dodge(.9),
                width = .5) +
  facet_wrap(~group, scales = "free") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_pattern_manual(values = c(Cells = "none", Medium = "stripe")) +
  labs(title = "Metabolites Levels (ng)")+
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = "white")))
ggsave("plots/010_Metabolites_concentration.png", h = 1500, w = 3000, units = "px")

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