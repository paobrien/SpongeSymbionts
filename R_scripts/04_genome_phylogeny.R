## Plot phylogeny of sponge and non-sponge associated genomes

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)
library(ape)
library(ggtree)
library(viridis)

## import data ----

# trees
trees <- readRDS(
  "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/tree_files_formatted"
)

# metadata
# get filepath
fp <- "~/Documents/R/Metagenomics/Data/Enrichm/Metadata/"

# get file names
meta_files <- list.files(fp, pattern = "meta_")

# batch import files
metadata <- list()
for(i in seq_along(meta_files)) {
  metadata[[i]] <- read.table(paste0(fp, meta_files[i]), 
                              header = T, sep = "\t", strip.white = T)
}

# name list
names(metadata) <- c("endo", "micro", "nitro", "spiro", "thermo")

# make tip labels
for(i in seq_along(metadata)) {
  metadata[[i]] <- metadata[[i]] %>%
    mutate(tip_label = paste0( Assembly_ID, " (", Isolation_short, ")"))
}

# check trees - individual taxa
#spiro
ggtree(trees$spiro, size=0.2) %<+% metadata$spiro + 
  geom_tiplab(size=2, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "grey", "blue")) +
  geom_nodelab(size=2) + xlim(0, 1.5)

#nitro
ggtree(trees$nitro, size=0.2) %<+% metadata$nitro + 
  geom_tiplab(size=2, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "grey", "blue")) +
  geom_nodelab(size=2) + xlim(0, 1.5) + ylim(0, 63)

#endo
ggtree(trees$endo, size=0.2) %<+% metadata$endo + 
  geom_tiplab(size=2, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "grey", "blue")) +
  geom_nodelab(size=2) + xlim(0, 0.6) 

#micro
ggtree(trees$micro, size=0.2) %<+% metadata$micro + 
  geom_tiplab(size=2, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "grey", "blue")) +
  geom_nodelab(size=2) + xlim(0, 0.9) 

#thermo
ggtree(trees$thermo, size=0.2) %<+% metadata$thermo + 
  geom_tiplab(size=2, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "grey", "blue")) +
  geom_nodelab(size=2) + xlim(0, 1) 

## plot all symbiont groups

# combine trees
trees$all <- bind.tree(trees$endo, trees$micro)
trees$all <- bind.tree(trees$all, trees$nitro)
trees$all <- bind.tree(trees$all, trees$spiro)
trees$all <- bind.tree(trees$all, trees$thermo)

# add symbiont name to metadata
metadata$endo$Symbiont <- rep("Endozoicomonadaceae", nrow(metadata$endo))
metadata$micro$Symbiont <- rep("Microtrichaceae", nrow(metadata$micro))
metadata$nitro$Symbiont <- rep("Nitrosopumiliaceae", nrow(metadata$nitro))
metadata$spiro$Symbiont <- rep("Spirochaetaceae", nrow(metadata$spiro))
metadata$thermo$Symbiont <- rep("Thermoanaerobaculaceae", nrow(metadata$thermo))

# combine metadata
metadata$all_groups <- do.call("rbind", metadata)

# remove outgroups
metadata$all_groups <- metadata$all_groups[!(metadata$all_groups$Group == "Outgroup"),]
trees$all <- keep.tip(trees$all, metadata$all_groups$Genome_ID)

# plot
# get cols
my_cols <- viridis(n = 5, alpha = 0.9)

# colour names
colour_names <- c("Nitrosopumilaceae", 
                  "Microtrichaceae", 
                  "Spirochaetaceae", 
                  "Thermoanaerobaculaceae", 
                  "Endozoicomonadaceae")

# create named vector
my_cols <- setNames(my_cols, colour_names)

# get nodes numbers for clades
ggtree(trees$all, size=0.2, layout = "circular") %<+% metadata$all_groups + 
  geom_text(aes(label=node), hjust=-.3, size = 2)

# get nodes to highlight
highlight_legend <- data.frame(
  node = c(281, 225, 341, 388, 213),
  group = colour_names
  )

highlight_legend$group <- factor(
  highlight_legend$group, 
  levels = colour_names)

# plot with ggtree and add highlights
p <- ggtree(trees$all, size=0.2, layout = "circular") %<+% metadata$all_groups +
  geom_tiplab2(size=1.9, linetype='dashed', align = T,
               linesize=.05, offset = .01, aes(label = tip_label, colour = Group)) + 
  scale_colour_manual(values = c("black", "blue"), guide = "none") +
  geom_nodelab(size=0.0001) + 
  geom_hilight(data = highlight_legend, 
               aes(node=node, fill=group),
               alpha = 0.6) +
  scale_fill_manual(values = my_cols) +
  guides(fill = guide_legend(title = "Taxonomy")) +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic", size = 13),
        legend.title = element_text(size = 15))
p

# save  
ggsave("all_symbionts_tree_v3.svg", 
       device = "svg", 
       width = 30, 
       height = 30, 
       dpi = 300, 
       units = "cm")
  













