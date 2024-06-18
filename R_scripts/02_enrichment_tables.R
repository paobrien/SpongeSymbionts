# Format enrichment data to make summary table

setwd("~/Documents/R/Metagenomics") 
library(tidyverse)


## load data ----
# cazy

# file path
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/CAZY/"

# list files
cazy_files <- list.files(fp, pattern = "CAZY_1")

cazy <- list()
for(i in seq_along(cazy_files)) {
  cazy[[i]] <- read.table(paste0(fp, cazy_files[[i]]), header = T, sep = ",", strip.white = T)
}
names(cazy) <- c("endo", "micro", "nitro", "spiro", "thermo")
lapply(cazy, head)

# ko
# file path
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/KO/"

# list files
ko_files <- list.files(fp, pattern = "KO_1")

ko <- list()
for(i in seq_along(ko_files)) {
  ko[[i]] <- read.table(paste0(fp, ko_files[[i]]), header = T, sep = ",", strip.white = T)
}
names(ko) <- c("endo", "micro", "nitro", "spiro", "thermo")
lapply(ko, head)


# pfam
# file path
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/PFAM/"

# list files
pfam_files <- list.files(fp, pattern = "PFAM_1")

pfam <- list()
for(i in seq_along(pfam_files)) {
  pfam[[i]] <- read.table(paste0(fp, pfam_files[[i]]), header = T, sep = ",", strip.white = T)
}
names(pfam) <- c("endo", "micro", "nitro", "spiro", "thermo")
lapply(pfam, head)


## Get total enriched

# function to get enriched
pos_enriched <- function(df) {
 return(filter(df, group_1_true > group_2_true & corrected_pvalue < 0.05) %>% nrow()) 
}

neg_enriched <- function(df) {
  return(filter(df, group_2_true > group_1_true & corrected_pvalue < 0.05) %>% nrow()) 
}

# positively enriched
# cazy
for(i in seq_along(cazy)) {
  print(names(cazy[i]))
  print(pos_enriched(cazy[[i]]))
}

# ko
for(i in seq_along(ko)) {
  print(names(ko[i]))
  print(pos_enriched(ko[[i]]))
}

# pfam
for(i in seq_along(pfam)) {
  print(names(pfam[i]))
  print(pos_enriched(pfam[[i]]))
}


# negatively enriched
# cazy
for(i in seq_along(cazy)) {
  print(names(cazy[i]))
  print(neg_enriched(cazy[[i]]))
}

# ko
for(i in seq_along(ko)) {
  print(names(ko[i]))
  print(neg_enriched(ko[[i]]))
}

# pfam
for(i in seq_along(pfam)) {
  print(names(pfam[i]))
  print(neg_enriched(pfam[[i]]))
}








