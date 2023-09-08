## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to process the scale cells and to
## generate plots for Figure 2 from processed data files available in the GEO repository GSE224695.
## Included here is also code for Figure 2 supplements 1-2

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e6)
})

source("Aman_Saunders_2023_utils.R")

# set genotype colors
geno.colors <- c("wt" = "#0077BB",
                 "nkt" = "#DDAA33",
                 "hypoth" = "#BB5566",
                 "bon" = "#009988")

# Load data ---------------------------------------------------------------

# all wildtype cells
wt_cds <- readRDS("data/GSM7029635_zskin_wildtype_cds.RDS")

# dermis cells
derm_cds <- readRDS("data/wt_dermis_cds.RDS")

# epidermis cells
epi_cds <-  readRDS("data/wt_epidermis_cds.RDS")

# Pre-processing dermis data -------------------------------------------

# do not run, processed data is loaded above
# an example of our dataset preprocessing

# 1. subset dermal data
derm_cds <- wt_cds[,colData(wt_cds)$cell_group == "Dermis"]

# 2. add cell cycle scores to the metadata
cc_list <- read.csv("data/cell_cycle_genes.csv", stringsAsFactors = F)

derm_cds <- estimate_cc_scores(derm_cds, 
                               gene_list1 = cc_list %>% filter(signature == "G1_S") %>% pull(id),
                               gene_list2 = cc_list %>% filter(signature == "G2_M") %>% pull(id))


# 3. process and regress cell cycle effects 
derm_cds = estimate_size_factors(derm_cds) %>% 
  detect_genes(derm_cds) %>% 
  preprocess_cds(derm_cds, num_dim = 30) %>% 
  reduce_dimension(derm_cds, preprocess_method = "PCA", 
                   umap.min_dist = 0.2, umap.n_neighbors = 20L)

# regress cell cycle effects from the dim reduction
derm_cds <- align_cds(derm_cds, 
                        residual_model_formula_str = "~proliferation_index", verbose = T)

# reduce dims and cluster
derm_cds <- reduce_dimension(derm_cds, 
                             preprocess_method = "Aligned", 
                             max_components = 2,
                             umap.min_dist = 0.2, 
                             umap.n_neighbors = 20L)

derm_cds <- cluster_cells(derm_cds, resolution = 1e-3)


# Subset and pre-process epidermal cells ----------------------------------

# do not run
# an example of our dataset preprocessing

# 1. subset epidermal cells
epi_cds <- wt_cds[,colData(wt_cds)$cell_group %in% c("Epidermis")]

# 2. preprocess, reduce dims and cluster cells
epi_cds <- detect_genes(epi_cds) %>% 
  estimate_size_factors() %>% 
  preprocess_cds(epi_cds, num_dim = 30) %>% 
  reduce_dimension(epi_cds, preprocess_method = "PCA", 
          umap.min_dist = 0.2, umap.n_neighbors = 20L)

epi_cds <- cluster_cells(epi_cds, resolution = 1e-3)

# Make plots for main figures -----------------------------------------------

# Figure 2a ---------------------------------------------------------------

# subset scale forming cells from all dermis
derm_filt <- derm_cds[,partitions(derm_cds) == 1]


# overlapping runx2 and sp7 gene expression in scale cells
plot_two_genes(derm_filt,"runx2b", "sp7") +
  theme(legend.position = "none")

ggsave(filename = "data/fig2a_scale_sp7_runx2b_co-expr_umap.png", 
       bg = "transparent", dpi = 600, 
       width = 4, height = 3.5)

# Figure 2f ---------------------------------------------------------------

# pseudotime analysis

derm_filt <- detect_genes(derm_filt) %>% 
  cluster_cells()

# root cells at the SFCs when prompted
derm_filt <- learn_graph(derm_filt) %>% 
  order_cells()

# plot
plot_cells(derm_filt, color_cells_by = "pseudotime", 
           show_trajectory_graph = F, alpha = 0.75, 
           cell_size = 1, cell_stroke = 0) +
  theme_void()

ggsave(filename = "data/fig2f_scale_pseudotime_umap.png", 
       bg = "transparent", dpi = 600, 
       width = 4, height = 3.5)

# Figure 2g ---------------------------------------------------------------

##  EMT signature analysis ##

# assign genes based on expression analysis and prior research 
epi_genes = data.frame("name" = c("kremen1", "cdh11", "cobl", "col11a1a"))
epi_genes$group = "epi"


emt_genes = data.frame(name = c("zeb2a", "zeb2b", "hdac4", "loxl2b", "snai2", "pdgfra"), 
                          group = "nc_mig")

sig_list = rbind(epi_genes, emt_genes)

sig_list <- sig_list %>% select(name, group) %>% 
  left_join(rowData(derm_filt) %>% 
              as.data.frame() %>% 
              select(id, gene_short_name), by = c("name" = "gene_short_name")) %>% 
  select("gene_id" = id, group)

# aggregate data to make scores across genes 
emt_mat = aggregate_gene_expression(cds = derm_filt, 
                                    gene_group_df = sig_list, 
                                    norm_method = "size_only",
                                    scale_agg_values = F)
sig_df = as.data.frame(t(emt_mat))
sig_df$Cell = rownames(sig_df)


# add migration ratio to metadata for plotting 
coldata_df <- colData(derm_filt) %>% 
  as.data.frame() %>% 
  left_join(sig_df, by = "Cell") %>% 
  dplyr::mutate(epi_mig_ratio = log10((nc_mig + 0.01)/(epi + 0.01)))

# scale and add bounds
coldata_df$epi_mig_scale <- scale(coldata_df$epi_mig_ratio)
coldata_df$epi_mig_scale[coldata_df$epi_mig_scale > 2] = 2
coldata_df$epi_mig_scale[coldata_df$epi_mig_scale < -2] = -2

# plot
ggplot(coldata_df, aes(umap1, umap2, color = epi_mig_scale)) +
  geom_point(size = 1, stroke = 0, alpha = 0.7, na.rm = T) +
  scale_color_distiller(palette = "RdBu", na.value = "grey80") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = "plots/fig3g_scale_emt-ratio_sig_expr_umap.png", 
       bg = "transparent", dpi = 600, 
       width = 4, height = 3.5)

# Figure 2i ---------------------------------------------------------------

## epidermal cell sub-analysis

# plot tp63 and krt4
plot_two_genes(epi_cds, "tp63", "krt4")

ggsave(filename = "plots/fig2i_epidermis_tp63-krt4_expr_umap_wl.png", 
       bg = "transparent", dpi = 750, width = 2, height = 2)


# Figure 2j ---------------------------------------------------------------

## stacked barplot of cells expressing two genes ##

# make a stacked barplot of the tp63/krt4 thresholding
epi_expr_df <- get_two_genes_df(epi_cds, "tp63", "krt4")

ggplot(epi_expr_df, aes(fill=expr, x=factor(cell_type_wt, levels = c("Basal Cell", "Suprabasal Cell", "Periderm")))) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = c("#4895ef", "#ef233c", "#9163cb", "grey90")) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))

ggsave(filename = "plots/fig2i_epidermis_tp63-krt4_stacked-barplot.pdf", 
       width = 3, height = 3.2)

# Figure supplements ------------------------------------------------------

# Figure 2-figure supplement 1b-h -----------------------------------------

## plot gene expression in scale cells

plot_genes <- c("enam", "spp1", "itga5", "chad", 
                "sp7", "scpp1", "axin2", "lef1", 
                "col10a1b", "col10a1a", "jag1a",
                "itgb3b", "itga5")

for (i in plot_genes[1:2]){
monocle3::plot_cells(derm_filt, 
                     genes = i, 
                     show_trajectory_graph = F,
                     label_cell_groups = F, 
                     cell_size = 0.75, 
                     cell_stroke = 0) +
  scale_color_viridis_c() +
  theme_void() +
  theme(legend.position = "right", # "none" for no legend
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(filename = paste0("data/fig2-S2_SFC_", i, "-expr_umap.png"),
       dpi = 750,
       height = 3,
       width = 3.75,
       bg = "transparent")
}



# Figure 2-figure supplement 2a -------------------------------------------

# plot proliferation index

plot_cells(derm_filt, 
           color_cells_by = "proliferation_index", 
           show_trajectory_graph = F, alpha = 0.75, 
           cell_size = .75, cell_stroke = 0) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = "data/fig2-S2b_scale_prolif-index_umap.png", 
       bg = "transparent", dpi = 750, 
       width = 2.5, height = 2)
