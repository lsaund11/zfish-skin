## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code generate plots for Figure 3 from 
## processed data files available in the GEO repository GSE224695.

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)

})

source("Aman_Saunders_2023_utils.R")

# Load data ---------------------------------------------------------------

# dermis cells
derm_cds <- readRDS("data/wt_dermis_cds.RDS")

# epidermis cells
epi_cds <-  readRDS("data/wt_epidermis_cds.RDS")

# Main figure plots ----------------------------------------------------------

# Figure 3b ---------------------------------------------------------------

## SCPP dotplot heatmap ##

# subset relevant celltypes 
epi_sub <- epi_cds[,partitions(epi_cds) %in% c(3, 4)] # basal cells + EMP epi
derm_sub <- derm_cds[,partitions(derm_cds) == 1 & 
                       colData(derm_cds)$cell_type_sub_v2 != "Pre-osteoclast"] # SFC pops

# collapse a few populations to broader types for this analysis

# epidermis (this also corrects a few stray cells)
colData(epi_sub)$cell_type_sub <- colData(epi_sub) %>% 
  as.data.frame() %>% 
  mutate(ct_fix = dplyr::recode(cell_type_sub, 
                                "Suprabasal Cell" = "Basal Cell",
                                "Periderm" = "Basal Cell",
                                "Inter-scale basal cells" = "Basal Cell",
                                "Proliferating BC" = "Basal Cell",
                                "Epidermal margin" = "Basal Cell")) %>% 
  pull(ct_fix)

# dermis
colData(derm_sub)$cell_type_sub <- colData(derm_sub) %>% 
  as.data.frame() %>% 
  mutate(ct_fix = dplyr::recode(cell_type_sub_v2, 
                                "Scale Pocket (eda, hhip+)" = "Pre-SFC")) %>% 
  pull(ct_fix)

# combine data and select genes (from separate analysis)
scpp_cds <- combine_cds(list(epi_sub, derm_sub), cell_names_unique = T)
scpp_cds <- detect_genes(scpp_cds)

scpp_genes <- c("enam",
               "ambn",
               "scpp5",
               "scpp7",
               "scpp8",
               "scpp1",
               "spp1",
               "lum",
               "chad",
               "sparc",
               "ostn",
               "bglap",
               "col10a1a",
               "col10a1b",
               "col2a1a",
               "col2a1b",
               "dlx3b",
               "dlx4a",
               "msx2a",
               "runx2b",
               "sp7")

expr_scpp <- rowData(scpp_cds) %>% 
  as.data.frame() %>% 
  filter(gene_short_name %in% scpp_genes) %>% 
  filter(num_cells_expressed > 50) %>% 
  pull(gene_short_name)

unique(colData(scpp_cds)$cell_type_sub)

ct_levs = c("Ameloblast", "Basal Cell", "SFC Sub-margin", 
            "SFC Margin", "Scale Radii Cell", "Pre-SFC", "Focus SFC")

colData(scpp_cds)$cell_type_sub <- factor(colData(scpp_cds)$cell_type_sub,
                                          levels = rev(ct_levs))

# plot dotplot heatmap
plot_genes_by_group(scpp_cds, scpp_genes,
                    group_cells_by = "cell_type_sub",
                    ordering_type = "none",
                    max.size = 4, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right", 
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 8, color = "black")) +
  coord_flip() +
  xlab("")

# save as pdf
ggsave("plots/fig3b_wt_basal-cell_SFC_scpp-gene-heatmap.pdf", 
       useDingbats = F, width = 5.8, height = 2.2)


# Figure 3c ---------------------------------------------------------------

## gene expression plots in basal cells ##

# ambn
plot_cells(epi_sub, genes = "ambn", label_cell_groups = F, cell_size = 0.75, cell_stroke = 0) +
  scale_color_viridis_c() +
  theme_void() #+
#theme(legend.position = "none")

ggsave(filename = "plots/fig2c_basalCell_ambn-expr_umap.png",
       dpi = 750,
       height = 2,
       width = 2.2,
       bg = "transparent")

# enam
plot_cells(epi_sub, genes = "enam", label_cell_groups = F, cell_size = 0.75, cell_stroke = 0) +
  scale_color_viridis_c() +
  theme_void() #+
#theme(legend.position = "none")

ggsave(filename = "plots/fig2c_basalCell_enam-expr_umap.png",
       dpi = 750,
       height = 2,
       width = 2.2,
       bg = "transparent")



# Figure 3g ---------------------------------------------------------------

## Scale forming cell subtype umap ##

sfc_cols <- c("SFC Sub-margin" = "#EA9100", "Pre-SFC" = "#b79ced",
                  "SFC Margin" = "#A71B4B",
                  "Scale Radii Cell" = "#fcde9c",
                  "Focus SFC" = "#f4a3a8")

derm_sub <- append_umap_coordinates(derm_sub)

derm_coldat <- colData(derm_sub) %>% 
  as.data.frame()

# with celltype labels
plot_cells(derm_sub, 
           color_cells_by = "cell_type_sub", 
           label_groups_by_cluster = F, 
           group_label_size = 4,
           show_trajectory_graph = F) +
  scale_color_manual(values = sfc_cols)


# without celltype labels
ggplot(derm_coldat) +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 1.2,
             stroke = 0) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = cell_type_sub),
             stroke = 0,
             size = 0.9) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_color_manual(values = sfc_cols)

ggsave(filename = "plots/fig2g_SFC_subtype_umap_nolab.png",
       bg = "transparent",
       width = 2.5,
       height = 2 ,
       dpi = 750)


# Figure 3i ---------------------------------------------------------------

## gene expression plots in SFCs ##

plot_genes <- c("spp1", "enam", "chad", "col10a1a")

for (gene in plot_genes){
    plot_cells(derm_sub, genes = gene, label_cell_groups = F, cell_size = 0.75, cell_stroke = 0) +
      scale_color_viridis_c() +
      theme_void() +
      theme(legend.position = "none", 
            strip.background = element_blank(),
            strip.text.x = element_blank())
    
    ggsave(filename = paste0("plots/fig2d_SFC_", gene, "-expr_umap.png"),
           dpi = 750,
           height = 2,
           width = 2.5,
           bg = "transparent")
}



