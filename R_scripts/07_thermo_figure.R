# Plot thermoanaerobaculaceae enrichment figure

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
  "~/Documents/R/Metagenomics/Data/Enrichm/Metadata/meta_thermo.tsv", 
  sep = "\t", 
  strip.white = T, 
  header = T
)

data_list <- readRDS(
  "~/Documents/R/Metagenomics/Data/Enrichm/Enrichment_figures/enrichment_data_formatted"
)

# remove trailing decimals from pfam colnames
colnames(data_list$pfam$thermo) <- gsub("\\..*", "", colnames(data_list$pfam$thermo))

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
trees$thermo$node.label <- gsub("d__.*", "", trees$thermo$node.label)
trees$thermo$node.label <- gsub(":.*", "", trees$thermo$node.label)

## plot tree ----

# check tree
p <- ggtree(trees$thermo, size=0.2) %<+% metadata + 
  geom_tiplab(aes(label = tip_label, colour = Group), 
              size=2, align=TRUE, linetype='dashed', 
              linesize=.05, offset = .1, show.legend=F) + 
  scale_colour_manual(values = c("black", "dark grey", "blue")) +
  geom_nodelab(size=1.8, hjust = 0.05) # + xlim(0, 1.5) + ylim(0, 50)
p

# check ankyrin genes
elp <- data_list$pfam$thermo[,grep("PF12796|PF13637|PF13857|PF00023",
                                   colnames(data_list$pfam$thermo))] 

# plot metadata
add_meta_gg(p, elp, "ANK", 0.5)

# check sod genes
sod <- data_list$ko$thermo[,grep("K04564|K04565|K00518",
                                colnames(data_list$ko$thermo))] 
# plot metadata
add_meta_gg(p, sod, "SOD", 0.5)


# check restriction enzymes
res <- data_list$ko$thermo[,grep("K01156|K01153|K01154|K07454|K01155",
                                 colnames(data_list$ko$thermo))] 

# plot metadata
add_meta_gg(p, res, "R-M", 0.5)

# check fibronectin & cadherin 

# for figures (enriched)
fib <- data_list$pfam$thermo[,grep("PF00041|PF00028",
                                   colnames(data_list$pfam$thermo))]

fib_t <- sqrt(fib)

# plot metadata
add_meta_gg(p, fib, "Binding", 0.5)


# check secretion system 
ss <- data_list$pfam$thermo[,grep("PF00482|PF00437|PF00263|PF05157|PF03958",
                                 colnames(data_list$pfam$thermo))]

ss_t <- sqrt(ss)

# plot metadata
add_meta_gg(p, ss, "SS", 0.5)


# check cazy
cazy <- data_list$cazy$thermo[, c("GH33", "CE10", "AA3", "GT20")]

# plot metadata
add_meta_gg(p, cazy, "CAZy", 0.5)

# check sulfur
sulf <- data_list$pfam$thermo[,grep("PF00884|PF14863|PF14864|PF13469|PF00685",
                                    colnames(data_list$pfam$thermo))]

sulf_t <- sqrt(sulf)

# plot metadata
add_meta_gg(p, sulf_t, "Sulfur", 0.5) + vexpand(.12, -1) + hexpand(.01, -1)

# check other enriched genes
og_k <- data_list$ko$thermo[,c("K03119", "K00387", "K01998", "K01995", "K01996")]

# plot metadata
add_meta_gg(p, og_k, "Other", 0.5)


# plot multiple on one tree

# edit column names
colnames(elp) <- c("Ank (PF00023)",
                   "Ank_2 (PF12796)",
                   "Ank_4 (PF13637)",
                   "Ank_5 (PF13857)")

colnames(res) <- c("hsdR (K01153)",
                   "hsdS (K01154)",
                   "(K01155)",
                   "res (K01156)",
                   "(K07454)")

colnames(sulf_t) <- c("Alkyl_sulf_C (PF14864)", 
                      "Alkyl_sulf_dimr (PF14863)", 
                      "Sulfatase (PF00884)", 
                      "Sulfotransfer_1 (PF00685)", 
                      "Sulfotransfer_3 (PF13469)")

colnames(ss) <- c("Secretin (PF00263)",
                  "Secretin_N (PF03958)",
                  "T2SSE (PF00437)",
                  "MshEN (PF05157)",
                  "T2SSF (PF00482)")

# set heatmap colours
mycol <- c("black", "purple", "white")

# plot tree
p1 <- add_meta_gg(p, elp, "ANK repeats", 0.5) + vexpand(.12, -1) + hexpand(.01, -1) # ANK repeats
p2 <- p1 + new_scale_fill()
p3 <- add_meta_gg(p2, res, "R-M system", 0.71) # restriction enzymes
p4 <- p3 + new_scale_fill()
p5 <- add_meta_gg(p4, cazy, "CAZy", 0.91)
p6 <- p5 + new_scale_fill()
p7 <- add_meta_gg(p6, sulf_t, "Sulfur (sqrt)", 1.12)
p8 <- p7 + new_scale_fill()
p9 <- add_meta_gg(p8, ss, "Secretion Systems", 1.325)
p9

p9 <- p9 + 
  annotate(geom="text", x=1.51, y=36, size = 4, label="ANK") +
  annotate(geom="text", x=1.72, y=36, size = 4, label="R-M") +
  annotate(geom="text", x=1.92, y=36, size = 4, label="CAZy") +
  annotate(geom="text", x=2.13, y=36, size = 4, label="Sulfur") +
  annotate(geom="text", x=2.33, y=36, size = 4, label="SS")
p9

# save
ggsave("thermo_tree.svg", 
       device = "svg",
       width = 28, 
       height = 25, 
       units = "cm",
       dpi = 300)







