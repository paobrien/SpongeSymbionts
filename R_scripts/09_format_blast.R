## Format tables to show which symbionts have taurine transporters

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)


## check genomes with gene ----

## import files ----

# batch import supp files
# get file path
fp <- "~/Documents/R/Metagenomics/Bin_stats/ch5_supp_tables/"

# get file names
supp_files <- list.files(fp, pattern = "Table")

# supp tables
supp_df <- list()
for(i in seq_along(supp_files)) {
  supp_df[[i]] <- read.table(paste0(fp, supp_files[i]), 
                              header = T, sep = "\t", strip.white = T)
}
names(supp_df) <- c("endo", "micro", "nitro", "spiro", "thermo")

# batch import taurine files
fp <- "~/Documents/Bioinformatics/PhD/Metagenomics/blast_taurine/significant_results_e03-100/" 

# get file names
tau_files <- list.files(fp, pattern = "genomes_with_taurine_unique_matches_e03")

# tau results
tau_df <- list()
for (i in seq_along(tau_files)) {
  tau_df[[i]] <- read.table(paste0(fp, tau_files[i]), 
                       header = F, sep = ";", strip.white = T)
}
names(tau_df) <- names(supp_df)



## format data ----

# add column names
for (i in seq_along(tau_df)) {
  names(tau_df[[i]]) <- c("Gene", "Genome")
}

# group by genome
tau_df <- lapply(tau_df, function(x) {
  x <- x %>%
    group_by(Genome) %>%
    summarize(Gene = paste(Gene, collapse = ";"))
}) 

# edit genome name
for (i in seq_along(tau_df)) {
  tau_df[[i]]$Genome <- gsub(".fna", "", tau_df[[i]]$Genome)
  tau_df[[i]]$Genome <- gsub("_genomic", "", tau_df[[i]]$Genome)
}

# combine dfs
for (i in seq_along(supp_df)) {
  supp_df[[i]] <- left_join(supp_df[[i]], tau_df[[i]], join_by("MAG_id" == "Genome"))
}


## make table of genomes with blast stats ----

## import files ----

# batch import significant results files
fp <- "~/Documents/Bioinformatics/PhD/Metagenomics/blast_taurine/significant_results_e03-100/" 

# get file names
result_files <- list.files(fp, pattern = "tblastx_significant_results")

# significant results
significant_results <- list()
for (i in seq_along(result_files)) {
  significant_results[[i]] <- read.table(paste0(fp, result_files[i]), 
                                         header = F, sep = "\t", strip.white = T)
}
names(significant_results) <- names(supp_df)

# batch import genome/contig with taurine files
# get file names
genome_files <- list.files(fp, pattern = "genomes_with_taurine_e03")

# genomes with taurine
genomes <- list()
for (i in seq_along(genome_files)) {
  genomes[[i]] <- read.table(paste0(fp, genome_files[i]), 
                             header = F, 
                             sep = ";", 
                             strip.white = T)
}
names(genomes) <- names(supp_df)


## format data ----

# significant results

# edit column names
for (i in seq_along(significant_results)) {
  names(significant_results[[i]]) <- c("Query_sequence_ID",	
                                       "Subject_sequence_ID", 	
                                       "Percent_identity",	
                                       "Alignment_length",	
                                       "Number_of_mismatches",	
                                       "Number_of_gap_openings",	
                                       "Start_position_in_query_sequence",	
                                       "End_position_in_query_sequence",	
                                       "Start_position_in_subject_sequence",	
                                       "End_position_in_subject_sequence",	
                                       "E-value",	
                                       "Bit_score")
}

# edit query seq
for (i in seq_along(significant_results)) {
  significant_results[[i]]$Query_sequence_ID <- gsub(";", "", significant_results[[i]]$Query_sequence_ID)
}

# genomes

# split contig/genome column
for (i in seq_along(genomes)) {
  genomes[[i]] <- separate(genomes[[i]], 
                           V2, 
                           c("Subject_sequence_ID", "Genome_ID"),
                           sep = " ")
}

# edit genome column
for (i in seq_along(genomes)) {
  names(genomes[[i]])[1] <- "Query_sequence_ID"
}

# edit genome ID
for (i in seq_along(genomes)) {
  genomes[[i]]$Genome_ID <- gsub(".fna", "", genomes[[i]]$Genome_ID)
  genomes[[i]]$Genome_ID <- gsub("_genomic", "", genomes[[i]]$Genome_ID)
}

# supp data
for (i in seq_along(genomes)) {
  supp_df[[i]] <- supp_df[[i]] %>%
    select(Genome_ID, Assembly_ID)
  genomes[[i]] <- left_join(genomes[[i]], supp_df[[i]])
}

# combine dataframes
for (i in seq_along(significant_results)) {
  significant_results[[i]] <- left_join(significant_results[[i]],
                                        genomes[[i]]
  )
}

# remove rows with NAs (genome not included)
for (i in seq_along(significant_results)) {
  significant_results[[i]] <- significant_results[[i]] %>%
    drop_na(Assembly_ID)
}

# replace genome id with assembly id
for (i in seq_along(significant_results)) {
  significant_results[[i]] <- significant_results[[i]] %>%
    select(!Genome_ID) %>%
    rename(Genome_ID = Assembly_ID)
}


## save ----

for (i in seq_along(significant_results)) {
  write.table(x = significant_results[[i]], 
              file = paste0("significant_blast_results_", names(significant_results)[i], ".tsv"), 
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t")
}









