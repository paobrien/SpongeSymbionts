# Plot microtrichaecae enrichment figure

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
  "~/Documents/R/Metagenomics/Data/Enrichm/Metadata/meta_micro.tsv", 
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

# make tip labels
metadata <- metadata %>%
  mutate(tip_label = paste0( Assembly_ID, " (", Isolation_short, ")"))

# edit node labels
trees$micro$node.label <- gsub("d__.*", "", trees$micro$node.label)
trees$micro$node.label <- gsub(":.*", "", trees$micro$node.label)
trees$micro$node.label <- gsub("f__.*", "", trees$micro$node.label)

## plot tree ----

# check tree
p <- ggtree(trees$micro, size=0.2) %<+% metadata + 
  geom_tiplab(aes(label = tip_label, colour = Group), 
              size=2, align=TRUE, linetype='dashed', 
              linesize=.05, offset = .1, show.legend = F) + 
  scale_colour_manual(values = c("black", "dark grey", "blue")) +
  geom_nodelab(size=1.8, hjust = 0.05) #+ xlim(0, 2) + ylim(0, 50)
p

# heatmap colours
mycol <- c("black", "purple")

# check ankyrin genes
elp <- data_list$pfam$micro[,grep("PF12796|PF13637|PF13857|PF00023",
                                  colnames(data_list$pfam$micro))] 
# plot metadata
add_meta_gg(p, elp, "ANK", 0.5) + vexpand(.05, -1)

# check restiction enzymes
res_k <- data_list$ko$micro[,grep("K03427|K01153|K19147|K01154|K01156",
                                  colnames(data_list$ko$micro))] 

# transform
res_k_t <- sqrt(res_k)

# heatmap colours
mycol <- c("black", "purple", "white")

# plot metadata
add_meta_gg(p, res_k, "R-M system", 0.5)


# check cas enzymes (pfam)
cas_p <- data_list$pfam$micro[,grep("PF01867|PF09617|PF09827|PF09609",
                                    colnames(data_list$pfam$micro))] 

# plot metadata
mycol <- c("black", "purple")
add_meta_gg(p, cas_p, "cas", 0.5)

# check cas enzymes (KO)
cas_k <- data_list$ko$micro[,grep("K15342|K07012|K19131|K09951|K19132",
                                  colnames(data_list$ko$micro))] 

# plot metadata
mycol <- c("black", "purple")
add_meta_gg(p, cas_k, "cas", 0.5)


# check host attachment: fibronectin & cadherin 
fib <- data_list$pfam$micro[,grep("PF00041|PF14310|PF05833|PF12733|PF00028",
                                  colnames(data_list$pfam$micro))]

fib_t <- sqrt(fib)

# heatmap colours
mycol <- c("black", "purple", "white")

# plot metadata
add_meta_gg(p, fib_t, "Binding (sqrt)", 0.5)


# check sulfatases/tranferases
sul <- data_list$pfam$micro[,grep("PF00884|PF03781|PF12411|PF13469|PF00685",
                                  colnames(data_list$pfam$micro))]

sul_t <- sqrt(sul)

# heatmap colours
mycol <- c("black", "purple", "white")

# plot metadata
add_meta_gg(p, sul_t, "Sulfur (sqrt)", 0.5)


# check additional genes
ag_ko <-  data_list$ko$micro[,grep("K01999|K01996|K01995|K01997|K01998",
                                   colnames(data_list$ko$micro))]

# plot metadata
add_meta_gg(p, ag_ko, "BCAA transport", 0.5)


## plot multiple on one tree ----

# edit column names
colnames(elp) <- c("Ank (PF00023)",
                   "Ank_2 (PF12796)",
                   "Ank_4 (PF13637)",
                   "Ank_5 (PF13857)")

colnames(res_k) <- c("hsdR (K01153)",
                      "hsdS (K01154)",
                      "res (K01156)",
                      "hsdM (K03427)", 
                      "mcrC (K19147)")

colnames(sul_t) <- c("Choline_sulf_C (PF12411)",
                     "FGE-sulfatase (PF03781)",
                     "Sulfatase (PF00884)",
                     "Sulfotransfer_1 (PF00685)",
                     "Sulfotransfer_3 (PF13469)")

colnames(fib_t) <- c("Cadherin (PF00028)", 
                     "Cadherin-like (PF12733)", 
                     "NFACT_N (PF05833)", 
                     "fn3 (PF00041)", 
                     "Fn3-like (PF14310)")

colnames(ag_ko) <- c("livG (K01995)",
                     "livF (K01996)",
                     "livH (K01997)",
                     "livM (K01998)",
                     "livK (K01999)")

# plot tree
# ankyrin repeats (set manual colours)
mycol <- c("black", "purple")
p1 <- gheatmap(p, elp, low = "black", high = "white", 
               color = "black", 
               offset = 0.5,
               colnames = T, 
               colnames_angle = 90,
               colnames_offset_y = 0.1,
               hjust=1,
               font.size=2.5,
               width = 0.2) + scale_fill_gradientn(name = "Ankyrin",
                                                   colours = mycol)  + 
  vexpand(.12, -1) + hexpand(.01, -1)

#p1 <- p1  + theme(legend.position = c(0.99, 1),
#                  legend.justification = c(0,1),
#                  legend.title=element_text(size = 12),legend.text=element_text(size=10))

# use function for remaining
mycol <- c("black", "purple", "white")
p2 <- p1 + new_scale_fill()
p3 <- add_meta_gg(p2, res_k, "R-M system", 0.675)
p4 <- p3 + new_scale_fill()
p5 <- add_meta_gg(p4, sul_t, "Sulfur (sqrt)", 0.85)
p6 <- p5 + new_scale_fill()
p7 <- add_meta_gg(p6, fib_t, "Binding (Sqrt)", 1.025)
p8 <- p7 + new_scale_fill()
p9 <- add_meta_gg(p8, ag_ko, "BCAA transport", 1.2)
p9

p9 <- p9 + 
  annotate(geom="text", x=1.35, y=60.25, size = 4, label="ANK") +
  annotate(geom="text", x=1.52, y=60.25, size = 4, label="R-M") +
  annotate(geom="text", x=1.695, y=60.25, size = 4, label="Sulfur") +
  annotate(geom="text", x=1.872, y=60.25, size = 4, label="Binding") +
  annotate(geom="text", x=2.04, y=60.25, size = 4, label="BCAA")
p9

# save
ggsave("micro_tree.svg", 
       device = "svg",
       width = 28, 
       height = 25, 
       units = "cm",
       dpi = 300)

