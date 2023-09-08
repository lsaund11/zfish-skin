## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to preprocess snRNA-seq data and to
## generate plots for Figure 1 from processed data files available in the GEO repository GSE224695.
## Included here is also code for Figure 1 supplements 2-4

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
geno.colors = c("wt" = "#0077BB",
                "nkt" = "#DDAA33",
                "hypoth" = "#BB5566",
                "bon" = "#009988")

# Load data --------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")
wt_cds <- readRDS("data/GSM7029635_zskin_wildtype_cds.RDS")

# Pre-processing for all genotype data -----------------------------------------------------

# do not run
# an example of our dataset preprocessing

all_cds <- detect_genes(all_cds) %>% 
  estimate_size_factors() %>% 
  preprocess_cds(num_dim = 50)

# use MNN alignment on genotype and background correct from metadata
all_cds <- align_cds(all_cds, alignment_group = "genotype", 
                 preprocess_method = "PCA",
                 residual_model_formula_str = "~ bg.loading.wt + bg.loading.nkt + 
                 bg.loading.bon + bg.loading.hypoth")

# UMAP dim reduction
all_cds <- reduce_dimension(all_cds, 
                       umap.min_dist = 0.1, 
                       umap.n_neighbors = 15L,
                       preprocess_method = "Aligned")

all_cds <- cluster_cells(all_cds, resolution = 1e-4)

# Pre-processing for wildtype data ----------------------------------------

## also don't run
## this has been done for the supplied dataset

# what are the genotypes called? 
unique(colData(all_cds)$genotype)

# subset data for wildtype
wt_cds <- all_cds[,colData(all_cds)$genotype == "wt"]

# continue with processing
wt_cds <- detect_genes(wt_cds) %>% 
  estimate_size_factors() %>% 
  preprocess_cds(num_dim = 50)

wt_cds <- align_cds(wt_cds, 
                    alignment_group = "sample")

wt_cds <- reduce_dimension(wt_cds, 
                           preprocess_method = "Aligned",
                           umap.min_dist = 0.15, 
                           umap.n_neighbors = 20L)

wt_cds <- cluster_cells(wt_cds, resolution = 2e-4)


# Manually assign cell types by cluster for wild type data -----------------

colData(wt_cds)$cluster <- clusters(wt_cds)

coldata_df <- colData(wt_cds) %>% 
  as.data.frame()

coldata_df <- coldata_df %>% 
  dplyr::mutate(cell_type_wt = case_when(
    cluster %in% c(2, 25) ~ "Suprabasal Cell",
    cluster %in% c(7, 6, 3) ~ "Periderm",
    cluster %in% c(10, 30) ~ "Muscle",
    cluster %in% c(5, 8, 16, 23) ~ "Basal Cell",
    cluster %in% c(1, 13, 14) ~ "Papillary Dermis",
    cluster == 28 ~ "Ionocyte",
    cluster == 32 ~ "NaR Ionocyte",
    cluster == 19 ~ "Xanthophore",
    cluster == 17 ~ "Iridophore",
    cluster == 24 ~ "Melanophore",
    cluster %in% c(15, 40, 31) ~ "Goblet Cell",
    cluster %in% c(4, 11, 41) ~ "SFC",
    cluster == 18 ~ "Hypodermis",
    cluster %in% c(12, 43, 36) ~ "Leukocyte",
    cluster == 9 ~ "Reticulate Dermis",
    cluster == 22 ~ "MLC",
    cluster == 34 ~ "PN/Schwann Cell",
    cluster == 27 ~ "Endothelial Cell",
    cluster == 20 ~ "PLL mantle cell",
    cluster == 26 ~ "PLL support cell",
    cluster == 29 ~ "PLL hair cell",
    TRUE ~ "Unknown"))

# add labels back to the cell data set 
colData(wt_cds)$cell_type_broad = coldata_df$cell_type_wt

plot_cells(wt_cds,
           color_cells_by = "cell_type_broad",
           label_groups_by_cluster = F, 
           group_label_size = 4)

# Add broad cell groups ---------------------------------------------------

coldata_df <- colData(wt_cds) %>% as.data.frame()

coldata_df <- coldata_df %>% 
  mutate(cell_group = case_when(cell_type_wt %in% c("Basal Cell", "Suprabasal Cell", "Periderm") ~ "Epidermis",
                                cell_type_wt %in% c("Reticulate Dermis", "SFC", "Papillary Dermis", "Hypodermis") ~ "Dermis",
                                cell_type_wt %in% c("Leukocyte", "MLC") ~ "Immune",
                                cell_type_wt %in% c("Iridophore", "Melanophore", "Xanthophore") ~ "Pigment",
                                cell_type_wt %in% c("PN/Schwann Cell") ~ "Neuron/Glia",
                                cell_type_wt %in% c("PLL mantle cell", "PLL support cell", "PLL hair cell") ~ "Lateral Line",
                                cell_type_wt == "Muscle" ~ "Muscle",
                                TRUE ~ "Other"))

colData(wt_cds)$cell_group <- coldata_df$cell_group

# Make plots for main figures -----------------------------------------------

# Define colors for cell types -------------------------------------------

# coloring groups for fig 1
celltype_cols_grey <- c("SFC" = "#A71B4B",
                       "Papillary Dermis" = "#EA9100",
                       "Reticulate Dermis" = "#F7E7A4",
                       "Epidermal margin" = "#B792BF",
                       "Hypodermis" = "#DECAA6",
                       "Melanophore" = "grey25",
                       "Xanthophore" = "#f2d338e5",
                       "Iridophore" = "#55A6CC",
                       "Periderm" = "#D4EBF9",
                       "Basal Cell" = "#D4EBF9",
                       "Suprabasal Cell" = "#D4EBF9",
                       "Unknown" = "grey80",
                       "Leukocyte" = "grey80",
                       #"Immune Cell" = "#007330",
                       "Goblet Cell" = "grey80",
                       "Muscle" = "#FFCCCC",
                       "PLL mantle cell" = "grey80", 
                       "PLL support cell" = "grey80", 
                       "PLL hair cell" = "grey80",
                       #"Hair Cell/Merkel Cell" = "#75ab7ee5",
                       "PN/Schwann Cell" = "grey80",
                       "Endothelial Cell" = "grey80", #D46D78
                       "MLC" = "grey80", #"#"
                       "Ionocyte" = "grey80",
                       "NaR Ionocyte" = "grey80")

# color all wildtype celltypes for supp figs
celltype_cols <- c("SFC" = "#A71B4B",
                  "Papillary Dermis" = "#EA9100",
                  "Reticulate Dermis" = "#F7E7A4",
                  "Hypodermis" = "#DECAA6",
                  "Epidermal margin" = "#B792BF",
                  "Melanophore" = "gray25",
                  "Xanthophore" = "#f2d338e5",
                  "Iridophore" = "#55A6CC",
                  "Periderm" = "#C1DAC5",
                  "Basal Cell" = "#007FA7",
                  "Suprabasal Cell" = "#17A690",
                  "Unknown" = "grey80",
                  "Leukocyte" = "#C4DE7B",
                  #"Immune Cell" = "#007330",
                  "Goblet Cell" = "#79B546",
                  "Muscle" = "#F3CFC6",
                  "PLL mantle cell" = "#1D6996", 
                  "PLL support cell" = "#1D6996", 
                  "PLL hair cell" = "#1D6996",
                  "PN/Schwann Cell" = "#994E95",
                  "Endothelial Cell" = "#D79883",
                  "MLC" = "#C796E2", #"#"
                  "Ionocyte" = "#584B9F",
                  "NaR Ionocyte" = "#584B9F")

# celltype colors for all genotypes (have a couple new ones)
all_celltype_cols <- c("SFC" = "#A71B4B",
                      "pre-SFC" = "#EA9100",
                      "Dermal mesenchyme" = "#F7E7A4",
                      "Myofibroblast" = "#F1D357",
                      "Hypodermis" = "#E5D3B5",
                      "Melanophore" = "gray25",
                      "Xanthophore" = "#f2d338e5",
                      "Iridophore" = "#55A6CC",
                      "Periderm" = "#C1DAC5",
                      "Basal Cell" = "#007FA7",
                      "Suprabasal Cell" = "#17A690",
                      "Unknown" = "grey80",
                      "Leukocyte" = "#C4DE7B",
                      "Leukocyte 2" = "#C4DE7B",
                      "Goblet Cell" = "#79B546",
                      "Goblet Cell 2" = "#79B546",
                      "Goblet Cell 3" = "#79B546",
                      "Muscle" = "#FFCCCC",
                      "Muscle 2" = "#FFCCCC",
                      "PLL mantle cell" = "#1D6996", 
                      "PLL support cell" = "#1D6996", 
                      "PLL hair cell" = "#1D6996",
                      "PN/Schwann Cell" = "#994E95",
                      "Schwann Cell" = "#994E95",
                      "Neuron" = "#994E95",
                      "Endothelial Cell" = "#D79883", 
                      "MLC" = "#C796E2",
                      "Ionocyte" = "#584B9F",
                      "NCC Ionocyte" = "#584B9F",
                      "Ionocyte (slc9a3.2+)" = "#584B9F",
                      "NCC Ionocyte" = "#584B9F",
                      "Pancreas" = "#C78006")

# Figure 1c ---------------------------------------------------------------

wt_coldata <- colData(wt_cds) %>% 
  as.data.frame()

# add the epidermal margin sub type annotation
em_cells <- wt_coldata %>% 
  filter(cell_type_sub == "Epidermal margin") %>% 
  pull(Cell)

wt_coldata[em_cells,]$cell_type_broad <- "Epidermal margin"

# plot
ggplot(wt_coldata %>%
         filter(!is.na(cell_type_broad)) %>%
         sample_n(size = dim(wt_coldata)[1])) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.75,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = cell_type_broad),
             stroke = 0,
             size = .6) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = celltype_cols_grey)

ggsave(filename = "plots/fig1c_wt_group-cols_umap.png",
       dpi = 750,
       height = 3,
       width = 3.8,
       bg = "transparent")

# Figure supplements ------------------------------------------------------

# Figure 1-figure supplement 2a-f -----------------------------------------

all_metadata <- colData(all_cds) %>% 
  as.data.frame() %>%
  mutate(genotype = factor(genotype, levels = c("wt", "nkt", "hypoth", "bon")))

## (a) UMIs per cell boxplot

ggplot(all_metadata %>% 
         filter(genotype == "wt") %>% 
         filter(!is.na(cell_type_broad)), 
       aes(x = cell_type_broad, y = log10(n.umi), fill = cell_type_broad)) +
  geom_boxplot(outlier.size = 0.5, size = 0.5, color = "black") +
  ylim(2, 4) +
  scale_fill_manual(values = all_celltype_cols) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()

ggsave(filename = "plots/fig1-s1a_all-geno_per-celltype_umi_boxplot.pdf", 
       width = 4, height = 3)

## (b) Genes per cell boxplot

ggplot(all_metadata %>% 
         filter(genotype == "wt") %>% 
         filter(!is.na(cell_type_broad)), 
       aes(x = cell_type_broad, y = num_genes_expressed, fill = cell_type_broad)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = all_celltype_cols) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()

ggsave(filename = "plots/fig1-s1b_all-geno_per-celltype_geneNum_boxplot.pdf", 
       width = 4, height = 3)

## (c) cell type umap

ggplot(all_metadata) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.7,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = cell_type_broad),
             stroke = 0,
             size = 0.5) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = all_celltype_cols)

ggsave(filename = "plots/fig1-s1c_all-geno_celltype_umap.pdf",
       dpi = 750,
       height = 3,
       width = 3,
       bg = "transparent")

# with cell type labels
plot_cells(all_cds, color_cells_by = "cell_type_broad", 
           label_groups_by_cluster = F,
           group_label_size = 3) +
  scale_color_manual(values = all_celltype_cols)

ggsave(filename = "plots/fig1-s1c_all-geno_celltype_umap_labs.pdf",
       dpi = 750,
       height = 3,
       width = 3,
       bg = "transparent")

## (d) genotype umap

ggplot(all_metadata) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.5,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = genotype),
             stroke = 0,
             size = .4) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = geno.colors)

ggsave(filename = "plots/fig1-s1c_all-geno_genotype_umap.pdf",
       dpi = 750,
       height = 3,
       width = 3,
       bg = "transparent")

## (e) genotype cell numbers

all_metadata %>% 
  select(Cell, genotype) %>% 
  ggplot(aes(x = genotype, fill = genotype)) +
  geom_bar(color = "black", size = 0.2) +
  scale_fill_manual(values = geno.colors) +
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  ggplot2::scale_y_continuous() +
  monocle3:::monocle_theme_opts()

ggsave(filename = "plots/fig1-s1e_all-geno_cell-counts_barplot.pdf", 
       width = 2, height = 3)

## (3) and genotype umi counts

ggplot(all_metadata, 
       aes(x = genotype, y = log10(n.umi), fill = genotype)) +
  geom_boxplot(outlier.size = 0.5, color = "black") +
  ylim(2, 4) +
  scale_fill_manual(values = geno.colors) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()

ggsave(filename = "plots/fig1-s1e_all-geno_umicounts_boxplot.pdf", 
       width = 2, height = 3)

## (f) genotype cell type compositions stacked bar plots

plot_df <- all_metadata %>% 
  group_by(genotype, cell_type_broad) %>% 
  tally(name = "count") %>% 
  ungroup() %>% 
  group_by(genotype) %>% 
  mutate(percent_per_cluster = count/sum(count) * 100) %>% 
  arrange(genotype)

ggplot(plot_df, aes(y = percent_per_cluster, x = genotype, fill = cell_type_broad)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) + 
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title.x=element_blank(), axis.text = element_text(color="Black"),
        axis.ticks = element_line(color="Black")) +
  scale_fill_manual(values = all_celltype_cols) +
  ylab("% Celltype") +
  monocle3:::monocle_theme_opts()

ggsave(filename = "plots/fig1-s1f_all-geno_per-celltype_stacked_barplot.pdf", 
         width = 2, height = 3, useDingbats = FALSE)


# Figure 1—figure supplement 3a ---------------------------------------------------------------

## umap colored by all cell types ##

ggplot(wt_coldata %>%
         filter(!is.na(cell_type_broad)) %>%
         sample_n(size = dim(wt_coldata)[1])) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.75,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = cell_type_broad),
             stroke = 0,
             size = .6) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = celltype_cols)

ggsave(filename = "plots/fig1-S3a_wt_broad-types_umap.png",
       dpi = 750,
       height = 3,
       width = 3.8,
       bg = "transparent")

# Figure 1—figure supplement 3b ---------------------------------------------------------------

## marker gene dotplot heatmap ##

# save genes and cell order
rev_genes <- c("tp63", "shha", "apoeb", "fgd5a", "muc5.1", "col6ab", "csf1b", 
              "kcnd2", "pnp4a", "ptprc", "tyrp1b", "adgrg11",
              "neb", "atp1b1b", "col6a3", "col2a1a", "col2a1b", "krt4",
              "sp7", "runx2b", "myo7aa", "mbpa", "pax7b", "cldna")

cell_order <- c("Basal Cell", "Epidermal margin", "Endothelial Cell", "Goblet Cell", 
               "Hypodermis", "Ionocyte", "Iridophore", "Leukocyte", "Melanophore", 
               "MLC", "Muscle", "NaR Ionocyte", "Papillary Dermis", 
               "Periderm", "SFC", "PLL hair cell", "PLL Mantle Cell", "PLL Support Cell", 
               "PN/Schwann Cell", "Xanthophore", "Reticulate Dermis", "Suprabasal Cell")

# clean and apply factors to order cells
plot_cds <- wt_cds[,colData(wt_cds)$cell_type_broad != "Unknown"]
colData(plot_cds)$cell_type_broad <- factor(colData(plot_cds)$cell_type_broad, levels = cell_order)

# plot
monocle3::plot_genes_by_group(plot_cds,
                    rev_genes,
                    ordering_type = 'none',
                    group_cells_by = "cell_type_broad",
                    max.size = 6, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right") +
  xlab("")

ggsave(filename = "plots/fig1-S3b_wt_dotplot_heatmap_wl.pdf", 
       useDingbats = F, width = 8, height = 6.5)

# Figure 1—figure supplement 4 ---------------------------------------------------------------

## Color wt cells by cluster from unsupervised leiden clustering ##

# plot clusters
monocle3::plot_cells(wt_cds, show_trajectory_graph = FALSE,
           color_cells_by = "cluster",
           label_cell_groups = T,
           group_label_size = 3,
           cell_size = 0.4) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = "plots/fig1-S4_zskin_wt_clusters_umap.pdf", 
       device = "pdf", width = 5, height = 4)


