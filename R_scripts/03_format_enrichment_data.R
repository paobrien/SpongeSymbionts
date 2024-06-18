## Format data for enrichment figures

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)
library(ape)
library(ggtree)


## Load and format files ----

## Load trees

# get filepath
fp <- "~/Documents/R/Metagenomics/Data/GTDB/"

# get file names
tree_files <- list.files(fp, pattern = "decorated.tree")

# batch import files
trees <- list()
for(i in seq_along(tree_files)) {
  trees[[i]] <- ape::read.tree(paste0(fp, tree_files[i]))
}
lapply(trees, function(x) length(x$tip.label))

# name list
names(trees) <- c("endo", "micro", "nitro", "spiro", "thermo")

# edit labels
for(i in seq_along(trees)) {
  trees[[i]]$tip.label <- gsub("RS_|GB_|_genomic", "", trees[[i]]$tip.label)
}

## load and format metadata

# get filepath
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/Metadata/"

# get file names
meta_files <- list.files(fp, pattern = "meta")

# batch import files
metadata <- list()
for(i in seq_along(meta_files)) {
  metadata[[i]] <- read.table(paste0(fp, meta_files[i]), 
                              header = T, sep = "\t", strip.white = T)
}
lapply(metadata, function(x) length(x$Genome_ID))

# name list
names(metadata) <- c("endo", "micro", "nitro", "spiro", "thermo")

# format files
for(i in seq_along(metadata)) {
  metadata[[i]]$Genome_ID <- gsub("_genomic", "", metadata[[i]]$Genome_ID)
  print(metadata[[i]]$Genome_ID)
  metadata[[i]]$Group <- gsub("1", "Sponge", metadata[[i]]$Group)
  metadata[[i]]$Group <- gsub("2", "Non-Sponge", metadata[[i]]$Group)
}


## Load data frequency tables

# cazy
# get filepath
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/CAZY/"

# get file names
cazy_files <- grep(list.files(fp), pattern = "allMags", invert=TRUE, value=TRUE) 

# batch import files
cazy <- list()
for(i in seq_along(cazy_files)) {
  cazy[[i]] <- read.table(paste0(fp, cazy_files[i]), 
                          header = T, sep = "\t", row.names = 1, 
                          check.names = F, strip.white = T)
}

# name list
names(cazy) <- c("endo", "micro", "nitro", "spiro", "thermo")


#ko
# get filepath
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/KO/"

# get file names
ko_files <- grep(list.files(fp), pattern = "_das_95", invert=TRUE, value=TRUE) 

# batch import files
ko <- list()
for(i in seq_along(ko_files)) {
  ko[[i]] <- read.table(paste0(fp, ko_files[i]), 
                        header = T, sep = "\t", row.names = 1, 
                        check.names = F, strip.white = T)
}

# name list
names(ko) <- c("endo", "micro", "nitro", "spiro", "thermo")


#pfam
# get filepath
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/PFAM/"

# get file names
pfam_files <- list.files(fp) 

# batch import files
pfam <- list()
for(i in seq_along(pfam_files)) {
  pfam[[i]] <- read.table(paste0(fp, pfam_files[i]), 
                          header = T, sep = "\t", row.names = 1, 
                          check.names = F, strip.white = T)
}

# name list
names(pfam) <- c("endo", "micro", "nitro", "spiro", "thermo")

# make data list
data_list <- list(cazy, ko, pfam)
names(data_list) <- c("cazy", "ko", "pfam")

## edit col names to match tip labels
# endo
for (i in seq_along(data_list)) {
  colnames(data_list[[i]]$endo) <- gsub("_genomic",
                                        "", colnames(data_list[[i]]$endo))
}

# micro
for (i in seq_along(data_list)) {
  colnames(data_list[[i]]$micro) <- gsub("_genomic",
                                         "", colnames(data_list[[i]]$micro))
}

# nitro
for (i in seq_along(data_list)) {
  colnames(data_list[[i]]$nitro) <- gsub("_genomic",
                                         "", colnames(data_list[[i]]$nitro))
}

# spiro
for (i in seq_along(data_list)) {
  colnames(data_list[[i]]$spiro) <- gsub("_20120800_E1X|_20120700_E2X|_20110800_E3D_genomic|_ASM381867v1|_genomic",
                                         "", colnames(data_list[[i]]$spiro))
}

# thermo
for (i in seq_along(data_list)) {
  colnames(data_list[[i]]$thermo) <- gsub("_genomic",
                                          "", colnames(data_list[[i]]$thermo))
}

## transpose data for heatmap
for(i in seq(data_list)) {
  data_list[[i]] <- lapply(data_list[[i]], t)
}

## Save ----

# Data for enrichment figures
saveRDS(trees, "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/tree_files_formatted")
saveRDS(metadata, "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/metadata_formatted")
saveRDS(data_list, "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/enrichment_data_formatted")


