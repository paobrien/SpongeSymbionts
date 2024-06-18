## Make supplementary tables of genome info 

setwd("~/Documents/R/Metagenomics") 

library(dplyr)

## Import data ----

## load and format metadata

group_names <- c("endo", "micro", "nitro", "spiro", "thermo")

# get filepath
fp <- "Data/Enrichm/Metadata/"

# get file names
meta_files <- list.files(fp, pattern = "meta")

# batch import files
metadata <- list()
for(i in seq_along(meta_files)) {
  metadata[[i]] <- read.table(paste0(fp, meta_files[i]), 
                              header = T, sep = "\t", strip.white = T)
}
lapply(metadata, function(x) length(x$MAG_id))

# name list
names(metadata) <- group_names

# format files
for(i in seq_along(metadata)) {
  metadata[[i]]$MAG_id <- gsub("_genomic", "", metadata[[i]]$MAG_id)
  print(metadata[[i]]$MAG_id)
  metadata[[i]]$Group <- gsub("1", "Sponge", metadata[[i]]$Group)
  metadata[[i]]$Group <- gsub("2", "Non-Sponge", metadata[[i]]$Group)
}

## load and format checkm data

# Import checkm results
filepath <- "Data/Checkm/"
checkm_files <- list.files(filepath, pattern = "long.tsv")[c(6:10)]
checkm <- lapply(checkm_files, 
                 function(x) read.table(paste0(filepath,x), 
                                        sep = "\t", 
                                        header = T, 
                                        strip.white = T))
# create list names
names(checkm) <- group_names


## load and format GTDB classification

# import gtdb tax
filepath <- "Data/Taxonomy/"
gtdb_files <- list.files(filepath, pattern = "_taxonomy")
gtdb <- lapply(gtdb_files,
               function(x) read.table(paste0(filepath,x), 
                                      sep = "\t", 
                                      header = T, 
                                      strip.white = T))

#names
names(gtdb) <- group_names

# edit mag id
for(i in seq_along(gtdb)) {
  gtdb[[i]]$MAG_id <- gsub("_genomic", "", gtdb[[i]]$MAG_id)
}

## load and format 16S classification

# import 16s silva tax
filepath <- "Data/Taxonomy/"

silva_files <- list.files(filepath, pattern = "_16s")
silva <- lapply(silva_files,
                function(x) read.csv(paste0(filepath,x), 
                                     sep = ",", 
                                     header = T, 
                                     strip.white = T))

#names
names(silva) <- group_names

# clean tax string and edit column and MAG id
for(i in seq_along(silva)) {
  silva[[i]]$Classification <- gsub("Root;_", "", silva[[i]]$Classification)
  silva[[i]]$MAG_id <- gsub("_genomic", "", silva[[i]]$MAG_id)
  names(silva[[i]])[2] <- "Classification_16S"
}

## Join and format data tables ----

# edit colnames to match metadata
for(i in seq_along(checkm)) {
  names(checkm[[i]])[1] <- "MAG_id"
}

# select columns to keep
for(i in seq_along(checkm)) {
  checkm[[i]] <- select(checkm[[i]], "MAG_id", "Completeness", 
                        "Contamination", "Genome_size_bp", "GC")
}

# edit genome size
for(i in seq_along(checkm)) {
  checkm[[i]]$Genome_size_bp <- checkm[[i]]$Genome_size_bp / 1000000
  colnames(checkm[[i]])[4] <- "Genome_size_Mbp"
}

# edit MAG ids
for(i in seq_along(checkm)) {
  checkm[[i]]$MAG_id <- gsub("_genomic", "", checkm[[i]]$MAG_id)
}

# join with metadata
df <- list()
for(i in seq_along(metadata)) {
  df[[i]] <- left_join(metadata[[i]], checkm[[i]], by = "MAG_id") 
  print(nrow(df[[i]]))
  print(nrow(metadata[[i]]))
}
names(df) <- c("endo", "micro", "nitro", "spiro", "thermo")

# left join gtdb tax
for(i in seq_along(df)) {
  df[[i]] <- left_join(df[[i]], gtdb[[i]], by = "MAG_id") 
  print(nrow(df[[i]]))
  print(nrow(gtdb[[i]]))
}

# left join 16S tax
for(i in seq_along(df)) {
  df[[i]] <- left_join(df[[i]], silva[[i]], by = "MAG_id") 
  print(nrow(df[[i]]))
  print(nrow(silva[[i]]))
}

## Save - export as tables ----

for(i in seq_along(df)) {
  write.table(df[[i]], 
              file = paste0("~/Documents/R/Metagenomics/Bin_stats/ch5_supp_tables/summary_", group_names[i], ".tsv"), 
              row.names = FALSE,
              col.names = TRUE, 
              sep = "\t")
}






