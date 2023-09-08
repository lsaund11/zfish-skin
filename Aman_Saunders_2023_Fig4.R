## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to generate plots for Figure 4 from 
## processed data files available in the GEO repository GSE224695.

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(scales)
})

source("Aman_Saunders_2023_utils.R")

# set genotype colors
geno.colors <- c("wt" = "#0077BB",
                "nkt" = "#DDAA33",
                "hypoth" = "#BB5566",
                "bon" = "#009988")

# Load data ---------------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")

# Figure 4a ---------------------------------------------------------------

## Eda and hypoth UMAPs ##

all_metadata <- colData(all_cds) %>% 
  as.data.frame()

# eda + wildtype
ggplot(all_metadata %>% 
         filter(genotype %in% c("wt", "nkt"))) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.7,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = genotype),
             stroke = 0,
             size = 0.5) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme_void() + 
  theme(legend.position = "right") +
  scale_color_manual(values = geno.colors)

ggsave(filename = "plots/fig4a_wt-eda_all-cells_umap.png",
       dpi = 750,
       height = 3,
       width = 3,
       bg = "transparent")

# hypothyroid + wildtype
ggplot(all_metadata %>% 
         filter(genotype %in% c("wt", "hypoth"))) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.7,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = genotype),
             stroke = 0,
             size = 0.5) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme_void() + 
  theme(legend.position = "right") +
  scale_color_manual(values = geno.colors)

ggsave(filename = "plots/fig4a_wt-hypoth_all-cells_umap.png",
       dpi = 750,
       height = 3,
       width = 3,
       bg = "transparent")


# Figure 4b  --------------------------------------------------------------

# Cell type abundance analysis ---------------------------------------------

counts_df <- all_metadata %>% 
  group_by(genotype, "cell_type" = cell_type_broad) %>% 
  filter(!is.na(cell_type)) %>% 
  dplyr::summarize(cells=n()) %>%
  ungroup() %>% 
  tidyr::pivot_wider(names_from = cell_type, 
                     values_from = cells, values_fill = c(0))

count_mat <- as.matrix(counts_df[,-1])
rownames(count_mat) <- counts_df$genotype

# grab meta data
meta_df <- all_metadata %>%
  select(genotype) %>% 
  distinct()
rownames(meta_df) <- meta_df$genotype

# make cell data set
count_mat <- count_mat[as.character(meta_df$genotype),] # same order

cell_cds <- new_cell_data_set(t(count_mat), 
                             cell_metadata=meta_df) %>% 
  preprocess_cds(num_dim = 10, 
                 norm_method="size_only", 
                 method = "PCA")

# extract sample size factors for normalization
sf_df <- size_factors(cell_cds) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column()
colnames(sf_df) = c("genotype" ,"size_factor")

# Assess celltype frequency across genotypes ------------------------------

# set factor levels
all_metadata$genotype <- factor(all_metadata$genotype, 
                               levels = c("wt", "nkt", "hypoth", "bon"))

# make normalized log-Fold change summary table across genotyoes
geno_counts <- all_metadata %>%
  dplyr::count(cell_type_broad, genotype, name = "count", .drop = FALSE) %>% 
  select("cell_type" = cell_type_broad, genotype, count) %>% 
  ungroup() %>%
  left_join(sf_df, by = "genotype") %>%
  mutate(norm_count = round(count/size_factor)) %>% 
  select(cell_type, genotype, norm_count) %>% 
  pivot_wider(names_from = genotype, 
              values_from = norm_count, 
              values_fill = list(norm_count = 0)) %>%
  mutate(nkt_log2fc = log2((nkt + 1)/wt)) %>% 
  mutate(bon_log2fc = log2((bon + 1)/wt)) %>% 
  mutate(hypoth_log2fc = log2((hypoth + 1)/wt)) %>% 
  select(cell_type, "Nkt" = nkt_log2fc, "Bon" = bon_log2fc, "hypoTH" = hypoth_log2fc)


# Plot heatmap  -----------------------------------------------------------

# select celltypes and clean up df
ct_sel <- all_metadata %>% 
  filter(cell_group %in% c("Dermis", "Epidermis")) %>%
  pull(cell_type_broad) %>% 
  unique()

hypo_nkt_plot <- geno_counts %>% 
  select(cell_type, Nkt, hypoTH) %>% 
  filter(cell_type %in% ct_sel) %>% 
  pivot_longer(cols = c(Nkt, hypoTH))

# plot
ggplot(hypo_nkt_plot, 
       aes(x = name, y = cell_type, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#3234a9", mid = "white", 
                       high= "darkred", na.value="black", name="") +
  theme(text = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.75),
        plot.margin=unit(c(.5,.5,.5,2.5),"cm"), legend.key.size = unit(0.5, 'cm'),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 0.25)) +
  coord_equal() +
  monocle3:::monocle_theme_opts()

ggsave("plots/fig4b_hypoth-eda_dermis-pigment_log2fc_heatmap.pdf", 
       width = 4, height  = 5, units = "in")

