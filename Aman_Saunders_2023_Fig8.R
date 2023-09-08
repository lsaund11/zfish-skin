## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code generate plots for Figure 8 and 
## the associated supplemental figure from processed data files available 
## in the GEO repository GSE224695.

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(pheatmap)
})

source("Aman_Saunders_2023_utils.R")

# set genotype colors
geno.colors = c("wt" = "#0077BB",
                "bon" = "#adeacf")

# Load data ---------------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")
wt_cds <- readRDS("data/GSM7029635_zskin_wildtype_cds.RDS")

# Main text figure plots --------------------------------------------------

# Figure 8b ---------------------------------------------------------------

## plot wt and bonaparte (bnc2) dermis, hypodermis and pigment cells ##

# subset wt+bon cells
wt_bon_cds <- all_cds[,colData(all_cds)$genotype %in% c("wt", "bon")]
wt_bon_cds <- wt_bon_cds[,partitions(wt_bon_cds) %in% c(5, 2, 8, 12, 11)]

# set factor levels and add log umi counts
colData(wt_bon_cds)$genotype <- factor(colData(wt_bon_cds)$genotype, levels = c("wt", "bon"))
colData(wt_bon_cds)$log.n.umi <- log10(colData(wt_bon_cds)$n.umi)

# check out cell types between genotypes
plot_cells(wt_bon_cds, 
           color_cells_by = "cell_type_broad", 
           label_groups_by_cluster = F, 
           group_label_size = 4) +
  facet_wrap(~genotype)

# extract metadata
coldat <- colData(wt_bon_cds) %>% 
  as.data.frame()

# plot for fig
ggplot(coldat) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = .6,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = genotype),
             stroke = 0,
             size = 0.4) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = geno.colors)

ggsave("plots/fig8b_wt_bon_derm-subset_umap.png",
       dpi = 500, width = 1.5, 
       height = 1.8, bg = "transparent")


# Figure 8c ---------------------------------------------------------------

# Cell type abundance analysis is the same as figure 4b
# that notebook has code to generate the `geno_counts` table used here

# Plot heatmap  -----------------------------------------------------------

# subsetting genotypes for heatmap
celltypes <- colData(all_cds) %>% 
  as.data.frame() %>% 
  filter(cell_group %in% c("Dermis", "Pigment")) %>%
  filter(cell_type_corr != "Myofibroblast") %>% 
  pull(cell_type_corr) %>% 
  unique()

bon_plot <- geno_counts %>% 
  select(cell_type, Bon) %>% 
  filter(cell_type %in% celltypes)

ggplot(bon_plot, 
       aes(x = 1, y = factor(cell_type, levels = c("Hypodermis", "pre-SFC",
                                                   "Dermal mesenchyme", "SFC", 
                                                   "Melanophore", "Iridophore", 
                                                   "Xanthophore")), fill = Bon)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#3234a9", mid = "white", high= "darkred", na.value="black", name="") +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(), text = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.75),
        plot.margin=unit(c(.5,.5,.5,2.5),"cm"), legend.key.size = unit(0.5, 'cm'),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 0.25),
        axis.ticks.x = element_blank()) +
  coord_equal() +
  monocle3:::monocle_theme_opts()

  ggsave("plots/fig8c_bon_dermis-pigment_cell-abund_heatmap.pdf", 
         width = 4, height  = 5, units = "in")
  

# Figure 8f ---------------------------------------------------------------

## Plotting pigment cell trophic factors expressed by dermis/hypodermis ##

# subset wildtype data
hd_pc_cds <- wt_cds[,colData(wt_cds)$cell_type_sub %in% c("Hypodermis", "Iridophore",
                                                          "Melanophore", "Xanthophore", 
                                                          "Myofibroblast")]

  
# make a cutoff by percentage of cells expressing per group
tf_order <- c("tjp1a", "csf1a", "csf1b", "asip1", 
               "edn3b", "aqp3a", "bnc2", "jam3b", 
               "kitlga", "igsf11", "asip2b", 
               "gdf5a", "alkal1", "alkal2a", "alkal2b")
  
# plot dotplot heatmap
monocle3::plot_genes_by_group(hd_pc_cds, 
                      markers = tf_order, 
                      group_cells_by = "cell_type_sub", 
                      max.size = 8) +
    theme(axis.title.x = element_blank())

ggsave("plots/fig8f_wt_bon_pc-trophic-factor_dotplot_heatmap.pdf", 
       width = 4, height = 4.5)


# Supplemental figure plots -----------------------------------------------

# Figure supplement 1b-c --------------------------------------------------

# edn3b
plot_cells(all_cds[,colData(all_cds)$genotype == "wt"], 
           genes = c("edn3b"), cell_size = 0.5, 
           label_cell_groups = F) +
  scale_color_viridis_c()
  
ggsave("plots/all-wt-sub_edn3b_expr_umap.png", 
       bg = "transparent", dpi = 700, 
       width = 3, height = 3.2, units = "in")

# csf1b
plot_umap_gene_expr(all_cds[,colData(all_cds)$genotype == "wt"], 
                    marker = "csf1b", cell_size = 0.5, 
                    label_cell_groups = F) +
  scale_color_viridis_c()

ggsave("plots/all-wt-sub_csf1b_expr_umap.png", 
       bg = "transparent", dpi = 700, 
       width = 3, height = 3.2, units = "in")

# Figure supplement 1d ----------------------------------------------------

# load differential expression data from DEA analysis
# this code is in the DEA notebook

# subset data
bon_degs <- wt_bon_degs %>% 
  filter(cell_type %in% c("Xanthophore", "Iridophore", "Melanophore")) %>% 
  group_by(cell_type, up_in, .drop = FALSE) %>% 
  tally(name = "deg_count") %>% 
  arrange(-deg_count, cell_type) %>% 
  ungroup()

# set factor levels
lvls <- bon_degs %>% 
  filter(up_in == "Bon") %>% 
  arrange(-deg_count) %>% 
  pull(cell_type)

bon_degs <- bon_degs %>% 
  filter(cell_type %in% lvls) %>% 
  mutate(cell_type = factor(cell_type, levels = lvls))

bon_degs$geno <- "bon"

# plot heatmap
ggplot(bon_degs, aes(x = cell_type, y = deg_count, fill = up_in)) +
  geom_bar(data = subset(bon_degs, up_in == "WT"), aes(y = deg_count), stat = "identity") +
  geom_bar(data = subset(bon_degs, up_in == "Bon"), aes(y = -deg_count), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        legend.background = element_rect(fill = "transparent", colour = NA)) + 
  scale_fill_manual(values = geno.colors) +
  monocle3:::monocle_theme_opts()

ggsave("plots/fig8-S1d_wt-bon_pigment-subset_deg_barplot.pdf", 
       bg = "transparent", 
       units = "in", width = 2.5, height = 3)
  