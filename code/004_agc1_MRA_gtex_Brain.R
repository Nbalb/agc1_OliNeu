library(DESeq2)
library(corto)
library(babelgene)
library(tibble)
library(dplyr)

setwd("data")
brainexpmat <- readRDS("data/gtex_Brain_expmat.rds")

descmatchr <- readRDS("data/002_res.rds") |> 
  as.data.frame() |> 
  filter(!duplicated(symbol))

rld <- readRDS("data/002_dds_rlog_normalized_agc1v2.rds")
agc1expmat <- assay(rld)

### Translate mouse genes to human gene
## Find most expressed orthologous gene in gtex for each mouse gene and select it to be the one to which the mouse symbol translates to
test_ort <- orthologs(rownames(descmatchr), species = "Mus musculus", human = F, top = T)
hsexpmat <- agc1expmat |> 
  as_tibble(rownames = "ensembl") |> 
  inner_join(test_ort, by = "ensembl") |> 
  select(human_symbol, starts_with(c("kd", "wt"))) |> 
  filter(!duplicated(human_symbol)) |>
  na.omit() |> 
  column_to_rownames("human_symbol") |> 
  as.data.frame()

### Load centroids
centroids <- read.delim("data/tfgenes_2020_09_11.txt", header = F)
centroids <- centroids[,2]

### GTEx Hippocampus Regulon ----
#Subset GTEx Brain to create a smaller regulon 
annot <- readRDS("data/gtex_sample_annotation.rds")
hippo <- annot[annot$SMTSD == "Brain - Hippocampus",]
hippoid <- hippo$SAMPID
hippoexpmat<- brainexpmat[,which(hippoid %in% colnames(brainexpmat))]
dim(hippoexpmat) #24274   123

### Create the regulon
fname <- "data/GTEx-Hippocampus_regulon.rda"
if(!file.exists(fname)){
  Hipporegulon <- corto(hippoexpmat, centroids, p = 1e-10, nbootstraps = 1000, verbose = T, nthreads = 6)
  saveRDS(Hipporegulon, file="data/GTEx-Hippocampus_regulon.rda")
}else(Hipporegulon <- readRDS("data/GTEx-Hippocampus_regulon.rda"))

### MRA on agc1 dataset
fname = "data/agc1-GTEx-Hippocampus-mra.rda"
if(!file.exists(fname)){
  ctr <- hsexpmat[,4:6]
  trt <- hsexpmat[,1:3]
  hippomra <- 
    mra(
      trt,
      ctr,
      regulon = Hipporegulon,
      minsize = 15,
      nthreads = 6,
      nperm = 1000,
      verbose = T
    )
  save(hippomra, file = "data/agc1-GTEx-Hippocampus-mra.rda")
}else(load(fname))

png("plots/mra-GTEx-Hippocampus.png", height = 10000, width = 6000, res=600)
mraplot(mraobj = hippomra, mrs = 10)
dev.off()

# ### GTEx Spinal Cord Regulon ----
# #Subset GTEx Brain to create a smaller regulon (GTEx Brain obteained got from Fede)
# load("data/annot.rda")
# unique(annot$SMTSD)
# spinalcord = annot[annot$SMTSD == "Brain - Spinal cord (cervical c-1)",]
# spinalcordid = spinalcord$SAMPID
# spinalcordexpmat = brainexpmat[,which(spinalcordid %in% colnames(brainexpmat))]
# dim(spinalcordexpmat) #24274   91
# 
# ### Create the regulon
# fname = "data/GTEx-Spinal-Cord_regulon.rda"
# if(!file.exists(fname)){
#   library(corto)
#   SCregulon <- corto(spinalcordexpmat,centroids,p=1e-10,nbootstraps = 1000,verbose = T,nthreads = 6)
#   saveRDS(SCregulon,file="data/GTEx-Spinal-Cord_regulon.rda")
# }else(SCregulon = readRDS(fname))
# 
# ### MRA on agc1 dataset
# fname = "data/agc1-GTEx-Spinal-Cord-mra.rda"
# if(!file.exists(fname)){
#   load("data/agc1-expmat.rda")
#   ctr<-hsexpmat[,4:6]
#   trt<-hsexpmat[,1:3]
#   SCmra <- mra(trt,ctr,regulon=SCregulon,minsize=15,nthreads=6,nperm=1000,verbose=T)
#   save(SCmra, file="data/agc1-GTEx-Spinal-Cord-mra.rda")
# }else(load(fname))
# 
# png("plots/mra-GTEx-Spinal-Cord.png", height = 7000, width = 6000, res=600)
# mraplot(mraobj = SCmra, mrs =10)
# dev.off()



### GTEx Frontal Cortex Regulon ----
#Subset GTEx Brain to create a smaller regulon (GTEx Brain obtained got from Fede)
unique(annot$SMTSD)
frontcortx = annot[annot$SMTSD == "Brain - Frontal Cortex (BA9)",]
frontcortxid = frontcortx$SAMPID
frontcortxexpmat = brainexpmat[,which(frontcortxid %in% colnames(brainexpmat))]
dim(frontcortxexpmat) #24274   129

### Create the regulon
fname = "data/GTEx-Frontal-Cortex_regulon.rda"
if(!file.exists(fname)){
  FCregulon <- corto(frontcortxexpmat,centroids,p=1e-10,nbootstraps = 1000,verbose = T,nthreads = 6)
  saveRDS(FCregulon, file = "data/GTEx-Frontal-Cortex_regulon.rda")
}else(FCregulon = readRDS(fname))

### MRA on agc1 dataset
fname = "data/agc1-GTEx-Frontal-Cortex-mra.rda"
if(!file.exists(fname)){
  FCmra <-
    mra(
      trt,
      ctr,
      regulon = FCregulon,
      minsize = 15,
      nthreads = 6,
      nperm = 1000,
      verbose = T
    )
  save(FCmra, file = "data/agc1-GTEx-Frontal-Cortex-mra.rda")
}else(load(fname))

png("plots/mra-GTEx-Frontal-Cortex.png", height = 10000, width = 6000, res=600)
mraplot(mraobj = FCmra, mrs = 10)
dev.off()

# Plot most interesting regulators, dlx1 and smarcc2
png("plots/dlx1_smarcc2_FC.png", height = 3000, width = 6000, res=600)
mraplot(mraobj = FCmra, mrs = c("DLX1", "SMARCC2"), pthr = 0.05)
dev.off()

png("plots/dlx1_smarcc2_Hippo.png", height = 3000, width = 6000, res=600)
mraplot(mraobj = hippomra, mrs = c("DLX1", "SMARCC2"), pthr = 0.05)
dev.off()
