## Plot spirochaete enrichment figure

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)

## Load data ----

trees <- readRDS(
  "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/tree_files_formatted"
)

metadata <- read.table(
  "~/Documents/R/Metagenomics/Data/Enrichm/Metadata/meta_spiro.tsv", 
  sep = "\t", 
  strip.white = T, 
  header = T
)

data_list <- readRDS(
  "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/enrichment_data_formatted"
)

# function to add metadata to ggtree
add_meta_gg <- function(p, df, legend_title, offset) {
  p1 <- gheatmap(p, df, low = "black", high = "white", 
                 color = "black", 
                 offset = offset,
                 colnames = T, 
                 colnames_angle = 90,
                 colnames_offset_y = 0.1,
                 hjust=1,
                 font.size=2.5,
                 width = 0.2)  +
    scale_fill_gradientn(name = legend_title,
                         colours = mycol)
  return(p1)
}

# set heatmap colours
mycol <- c("black", "purple", "white")

# make tip labels
metadata <- metadata %>%
  mutate(tip_label = paste0( Assembly_ID, " (", Isolation_short, ")"))

# edit node labels
trees$spiro$node.label <- gsub("d__.*", "", trees$spiro$node.label)
trees$spiro$node.label <- gsub(":.*", "", trees$spiro$node.label)
trees$spiro$node.label <- gsub("f__.*", "", trees$spiro$node.label)
trees$spiro$node.label <- gsub("Alkalispirochaetaceae", "", trees$spiro$node.label)


## plot tree ----

# check tree
p <- ggtree(trees$spiro, size=0.2) %<+% metadata + 
  geom_tiplab(aes(label = tip_label, colour = Group), 
              size=2, align=TRUE, linetype='dashed', 
              linesize=.05, offset = .1, show.legend = F) + 
  scale_colour_manual(values = c("black", "dark grey", "blue")) +
  geom_nodelab(size=2, hjust = 0.11) 
p

# check cazy data 

# subset cazy genes to plot
cazy <- data_list$cazy$spiro[,c("GH29", "GH33", "GH51", "GT1","GT94")]

# transform data
cazy_t <- sqrt(cazy)

# add heatmap to tree
add_meta_gg(p, cazy_t, "CAZy (sqrt)", 0.5) + vexpand(.05, -1)


## check ankyrin repeats
ank <- data_list$pfam$spiro[,grep("PF13637|PF12796|PF00023|PF13857",
                                  colnames(data_list$pfam$spiro))] 

# plot metadata
add_meta_gg(p, ank, "ANK repeats", 0.5) + vexpand(.05, -1)

## check sulfur enzymes
slf <- data_list$pfam$spiro[,grep("PF14269|PF05935|PF00884|PF13469|PF03781",
                                  colnames(data_list$pfam$spiro))] 

# transform
slf_t <- sqrt(slf)

# plot metadata
add_meta_gg(p, slf_t, "Sulfur (sqrt)", 0.5) + vexpand(.05, -1)

# check restriction enzymes
rm_sys <- data_list$ko$spiro[, c("K07454", "K01156", "K01154", "K01153", "K03427")]

# plot metadata
add_meta_gg(p, rm_sys, "R-M system", 0.5) + vexpand(.05, -1)

## check ureases
urea <- data_list$ko$spiro[, c("K01428", "K01430", "K03188", "K01429", "K03189")]

# change heatmap colours (as only 0-1)
mycol <- c("black", "purple")

# plot metadata
add_meta_gg(p, urea, "Ureases", 0.5) + vexpand(.05, -1)

## plot multiple on one tree ----

# set heatmap colours
mycol <- c("black", "purple", "white")

# edit column names
colnames(ank) <- c("Ank (PF00023)",
                   "Ank_2 (PF12796)",
                   "Ank_4 (PF13637)",
                   "Ank_5 (PF13857)")

colnames(rm_sys) <- c("(K07454)",
                      "res (K01156)",
                      "hsdS (K01154)",
                      "hsdR (K01153)", 
                      "hsdM (K03427)")

colnames(cazy_t) # already gene codes
colnames(slf_t) <- c("Arylsulfotrans (PF05935)",
                     "Arylsulfotran_2 (PF14269)",
                     "FGE-sulfatase (PF03781)",
                     "Sulfatase (PF00884)",
                     "Sulfotransfer_3 (PF13469)")

colnames(urea) <- c("ureC (K01428)", 
                    "ureA (K01430)", 
                    "ureF (K03188)", 
                    "ureB (K01429)", 
                    "ureG (K03189)")

# plot combined tree
p1 <- add_meta_gg(p, ank, "ANK repeats", 0.5) + vexpand(.12, -1)
p2 <- p1 + new_scale_fill()
p3 <- add_meta_gg(p2, rm_sys, "R-M system", 0.75)
p4 <- p3 + new_scale_fill()
p5 <- add_meta_gg(p4, cazy_t, "CAZy (sqrt)", 1)
p6 <- p5 + new_scale_fill()
p7 <- add_meta_gg(p6, slf_t, "Sulfur (sqrt)", 1.25)
p8 <- p7 + new_scale_fill()
p9 <- gheatmap(p8, urea, low = "black", high = "white", 
         color = "black", 
         offset = 1.50,
         colnames = T, 
         colnames_angle = 90,
         colnames_offset_y = 0.1,
         hjust=1,
         font.size=2.5,
         width = 0.2)  +
  scale_fill_gradientn(name = "Ureases",
                       colours = c("black", "purple"))
p9 <- p9 + 
  annotate(geom="text", x=1.69, y=51.1, size = 4, label="ANK") +
  annotate(geom="text", x=1.94, y=51.1, size = 4, label="R-M") +
  annotate(geom="text", x=2.185, y=51.1, size = 4, label="CAZy") +
  annotate(geom="text", x=2.43, y=51.1, size = 4, label="Sulfur") +
  annotate(geom="text", x=2.685, y=51.1, size = 4, label="Ureases")
p9


# save
ggsave("spiro_tree.svg", 
       device = "svg",
       width = 28, 
       height = 25, 
       units = "cm",
       dpi = 300)
