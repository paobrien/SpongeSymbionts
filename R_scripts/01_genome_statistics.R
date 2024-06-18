## Get genome statistics and plot GC content and genome size

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)
library(reshape2)
library(scales)
library(forcats)
library(ggpubr)
library(viridis)

## Import data ----

# import checkm results
filepath <- "Data/Checkm/"
checkm_files <- list.files(filepath, pattern = "long.tsv")[c(6:10)]
checkm <- lapply(checkm_files, 
                 function(x) read.table(paste0(filepath,x), 
                                        sep = "\t", 
                                        header = T, 
                                        strip.white = T))
# create list names
names_sym <- c("endo", "micro", "nitro", "spiro", "thermo")
names(checkm) <- names_sym

# import group metadata
filepath <- "Data/Enrichm/Metadata/"
meta_files <- list.files(filepath, pattern = "meta")
metadata <- lapply(meta_files, 
                   function(x) read.table(paste0(filepath,x), 
                                          sep = "\t", 
                                          header = T, 
                                          strip.white = T))

# edit col name to match checkm file
for(i in seq_along(metadata)) {
  names(metadata[[i]])[1] <- "Bin_Id"
}

# create list names
names(metadata) <- names_sym

# edit genome names
#for (i in seq_along(checkm)) {
#  checkm[[i]]$Bin_Id <- gsub("_genomic", "",checkm[[i]]$Bin_Id)
#}

# remove columns not in analysis (those removed through dereplication)
for(i in seq_along(checkm)) {
  checkm[[i]] <- checkm[[i]][(checkm[[i]]$Bin_Id %in% metadata[[i]]$Bin_Id),]
} 

# check bins are the same in both dataframes
for(i in seq_along(metadata)) {
  print(paste(names(metadata[i]), 
              setdiff(metadata[[i]]$Bin_Id, checkm[[i]]$Bin_Id)))
}

# add grouping column and format
for(i in seq_along(checkm)) {
  checkm[[i]] <- left_join(checkm[[i]], metadata[[i]], by = "Bin_Id")
  checkm[[i]]$Group <- gsub("1", "Sponge", checkm[[i]]$Group)
  checkm[[i]]$Group <- gsub("2", "Non-Sponge", checkm[[i]]$Group)
  checkm[[i]]$Genome_size_Mbp <- checkm[[i]]$Genome_size_bp / 1000000
}
head(checkm$endo)


# get mean and SE of for each group
sum_stats <- list()
for(i in seq_along(checkm)) {
  sum_stats[[i]] <- group_by(checkm[[i]], Group) %>% 
    summarise(mean_size = mean(Genome_size_Mbp, na.rm = T), 
              stdev_size = sd(Genome_size_Mbp, na.rm = T), 
              sterr_size = sd(Genome_size_Mbp, na.rm = T)/sqrt(n()),
              min_size = min(Genome_size_Mbp, na.rm = T),
              max_size = max(Genome_size_Mbp, na.rm = T),
              mean_gc = mean(GC, na.rm = T),
              stdev_gc = sd(GC, na.rm = T),
              sterr_gc = sd(GC, na.rm = T)/sqrt(n()),
              min_gc = min(GC, na.rm = T),
              max_gc = max(GC, na.rm = T)) %>% as.data.frame
}

names(sum_stats) <- names_sym
sum_stats


## Save ----
saveRDS(checkm, "~/Documents/R/Metagenomics/Data/Genome_summary_stats/checkm_file_list")
saveRDS(sum_stats, "~/Documents/R/Metagenomics/Data/Genome_summary_stats/summary_stats")

checkm <- readRDS("~/Documents/R/Metagenomics/Data/Genome_summary_stats/checkm_file_list")
sum_stats <- readRDS("~/Documents/R/Metagenomics/Data/Genome_summary_stats/summary_stats")


## ttest genome size and gc content ----

# Look for statistic difference between groups

# check for normality
shapiro_size <- list()
shapiro_gc <- list()
for(i in seq_along(checkm)) {
  shapiro_size[[i]] <- shapiro.test(checkm[[i]]$Genome_size_Mbp)
  shapiro_gc[[i]] <- shapiro.test(checkm[[i]]$GC)
}
shapiro_size
shapiro_gc

# using welches ttest due to non-normality
# first remove outgroup
for(i in seq_along(checkm)) {
  checkm[[i]] <- checkm[[i]][!(checkm[[i]]$Group=="Outgroup"),] 
}

# run ttest
ttest_size <- list()
ttest_gc <- list()
for(i in seq_along(checkm)) {
  ttest_size[[i]] <- t.test(Genome_size_Mbp~Group, checkm[[i]], var.equal = FALSE)
  ttest_gc[[i]] <- t.test(GC~Group, checkm[[i]], var.equal = FALSE)
}

names(ttest_size) <- names_sym
names(ttest_gc) <- names_sym

# save
saveRDS(ttest_size, "~/Documents/R/Metagenomics/Data/Genome_summary_stats/ttest_size")
saveRDS(ttest_gc, "~/Documents/R/Metagenomics/Data/Genome_summary_stats/ttest_gc")

ttest_size <- readRDS("~/Documents/R/Metagenomics/Data/Genome_summary_stats/ttest_size")
ttest_gc <- readRDS("~/Documents/R/Metagenomics/Data/Genome_summary_stats/ttest_gc")


## Create table of summary stats ----
# drop outgroup row
for(i in seq_along(sum_stats)) {
  sum_stats[[i]] <- sum_stats[[i]][!(sum_stats[[i]]$Group=="Outgroup"),]
}

# combine lists
sum_stats <- do.call(rbind, sum_stats)

# add microbe tax column
sum_stats$Family <- c("Endo", "Endo",
                      "Micro", "Micro",
                      "Nitro", "Nitro",
                      "Spiro", "Spiro",
                      "Thermo", "Thermo")
# reorder
sum_stats <- select(sum_stats, Group, Family, everything())


# export table
write.table(sum_stats, 
            file = "~/Documents/R/Metagenomics/Data/Genome_summary_stats/summary_stats_all.tsv", 
            row.names = FALSE,
            col.names = TRUE, 
            sep = "\t")


## Plot gc vs. genome size (Mbp) ----

## note: load from Symbiont Save above

# remove outgroup
for (i in seq_along(checkm)) {
  checkm[[i]] <- checkm[[i]][!(checkm[[i]]$Group == "Outgroup"),]
}

# add symbiont column
for (i in seq_along(checkm)) {
  checkm[[i]]$Symbiont <- rep(names(checkm)[i], nrow(checkm[[i]]))
}

# subset columns
for (i in seq_along( checkm)) {
  checkm[[i]] <- select(checkm[[i]], Bin_Id, 
                        Completeness, 
                        Contamination, 
                        GC, 
                        Isolation_source, 
                        Genome_size_Mbp,
                        Group,
                        Coding_density,
                        Symbiont)
}

# combine dfs 
checkm$all_groups <- do.call("rbind", checkm)


# plot geonome size vs GC
p_endo <- ggplot(checkm$endo, aes(x = Genome_size_Mbp, y = GC, fill = Group)) + 
  geom_point(alpha = 0.8, size = 5, shape = 21) +
  scale_fill_viridis(discrete = T) +
  xlim(0, 8.5) +
  ylim(20, 75) +
  theme_bw()

p_micro <- ggplot(checkm$micro, aes(x = Genome_size_Mbp, y = GC, fill = Group)) + 
  geom_point(alpha = 0.8, size = 5, shape = 21) +
  scale_fill_viridis(discrete = T) +
  xlim(0, 8.5) +
  ylim(20, 75) +
  theme_bw()

p_nitro <- ggplot(checkm$nitro, aes(x = Genome_size_Mbp, y = GC, fill = Group)) + 
  geom_point(alpha = 0.8, size = 5, shape = 21) +
  scale_fill_viridis(discrete = T) +
  xlim(0, 8.5) +
  ylim(20, 75) +
  theme_bw()

p_spiro <- ggplot(checkm$spiro, aes(x = Genome_size_Mbp, y = GC, fill = Group)) + 
  geom_point(alpha = 0.8, size = 5, shape = 21) +
  scale_fill_viridis(discrete = T) +
  xlim(0, 8.5) +
  ylim(20, 75) +
  theme_bw()

p_thermo <- ggplot(checkm$thermo, aes(x = Genome_size_Mbp, y = GC, fill = Group)) + 
  geom_point(alpha = 0.8, size = 5, shape = 21) +
  scale_fill_viridis(discrete = T) +
  xlim(0, 8.5) +
  ylim(20, 75) +
  theme_bw()

# plot boxplot of genome size and GC

p_gs <- ggplot(checkm$all_groups, aes(x=Symbiont, y=Genome_size_Mbp, fill=Group)) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_viridis(discrete = T) +
  theme_bw() + 
  theme(legend.position="none")

p_gc <-  ggplot(checkm$all_groups, aes(x=Symbiont, y=GC, fill=Group)) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_viridis(discrete = T) + 
  theme_bw() +
  theme(legend.position="none")

p_all <- ggarrange(p_gs+theme(axis.text.x = element_blank(), 
                              axis.title.x = element_blank(),
                              axis.ticks.x = element_blank()), 
                   p_gc, nrow = 1, ncol = 2)

# boxplot of genome completeness
p_comp <- ggplot(checkm$all_groups, aes(x=Symbiont, y=Completeness, fill=Group)) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_viridis(discrete = T) +
  ylim(70, NA) +
  theme_bw() + 
  theme(legend.position="none")


# plot all on one figure
ggarrange(p_endo +
            ylab("") + ggtitle("Endozoicomonadaceae") +
            theme(axis.title.x = element_blank(), 
                  plot.title = element_text(face="italic"),
                  legend.text=element_text(size=11), legend.title=element_text(size=13)), 
          p_micro + 
            ylab("")  + ggtitle("Microtrichaceae") +
            theme(axis.title.x = element_blank(), plot.title = element_text(face="italic")), 
          p_nitro + 
            ggtitle("Nitrosopumiliaceae") +
            theme(axis.title.x = element_blank(), plot.title = element_text(face="italic")),
          p_spiro + 
            ylab("") + ggtitle("Spirochaetaceae") +
            theme(axis.title.x = element_blank(), plot.title = element_text(face="italic")),
          p_thermo + 
            ylab("") + ggtitle("Thermoanaerobaculaceae") +
            theme(plot.title = element_text(face="italic")),
          p_comp + ggtitle(""), 
          heights = c(0.9, 0.9, 1.0),
          nrow = 3, ncol = 2,  common.legend = T, legend = "right")

# save
ggsave("genome_plot_test2.pdf", device = "pdf", width = 20, height = 25, units = "cm", dpi = 300)











