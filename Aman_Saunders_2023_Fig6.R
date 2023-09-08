## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to generate plots for Figure 6 
## from processed data files available in the GEO repository GSE224695.

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
geno.colors <- c("wt" = "#0077BB",
                 "nkt" = "#DDAA33",
                 "hypoth" = "#BB5566",
                 "bon" = "#009988")

# set heatmap colors
hm_colors <-c("white", "#fad2d2", "#ce4257", "#720026")


# Load data ---------------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")
wt_cds <- readRDS("data/GSM7029635_zskin_wildtype_cds.RDS")

# Figure 6a ---------------------------------------------------------------

## HypoTH DEG heatmap ##

hypoth_dea_df <- fread(file = "data/wt-hypoth_all-celltype_fit-model_res.csv",
                       sep = ",", data.table = F, na.strings = "NA")

# select cell types
celltypes <- c("Basal Cell", "Periderm",
              "Suprabasal Cell",
              "Dermal mesenchyme", "Hypodermis")

# get list of DEGs for heatmap
deg_ids <- hypoth_dea_df %>%
  filter(cell_type %in% celltypes) %>% 
  filter(q_value < 0.01) %>%
  filter(normalized_effect >= abs(1)) %>% 
  pull(id) %>%
  unique()

length(deg_ids) # 836 degs

# take a look
hypoth_dea_df %>% filter(q_value < 0.01) %>%
  arrange(-normalized_effect) %>% 
  filter(cell_type == "Basal Cell") %>% 
  head()

# check distribution of test betas
hist(hypoth_dea_df$normalized_effect, xlim = c(-5, 5))

# make matrix of the test betas
beta_wide_df <- hypoth_dea_df %>%
  dplyr::filter(id %in% deg_ids) %>%
  dplyr::filter(cell_type %in% celltypes) %>%
  dplyr::distinct() %>%
  dplyr::select(id, normalized_effect, cell_type) %>% 
  tidyr::pivot_wider(
    names_from = cell_type,
    values_from = normalized_effect,
    values_fill = c(0))

beta_mat <- as.matrix(beta_wide_df[,-1])
rownames(beta_mat) <- beta_wide_df$id
beta_mat[beta_mat > 3] <- 3
beta_mat[beta_mat < -3] <- -3

# reorder column names and check
beta_mat <- beta_mat[,celltypes]
colnames(beta_mat)
beta_mat[1:5, 1:5]

BuRd <- c('#2166AC', '#4393C3', '#92C5DE', '#D1E5F0', '#F7F7F7', "#fad2d2", "#f08080", "#ce4257", "#720026")

# plot heatmap
pdf("plots/fig6a_hypoth_DEG_heatmap_across-celltypes.pdf", width = 5, height = 6, bg = "transparent", useDingbats = F)           
pheatmap(beta_mat,
         color = colorRampPalette(BuRd)(100),
         cluster_rows = T,
         cluster_cols = F, border_color = "NA",
         show_rownames = F, show_colnames = T, clustering_method = "ward.D2")
dev.off()

# Figure 6b ---------------------------------------------------------------

## Pathway ligand heatmap ##

# load DEG ligand resuts for all genotypes 
dea_ligand_res <- fread("data/pathway_ligand_DEG_res.csv", sep = ",", data.table = F)
  
# select celltypes
celltypes <- c("Basal Cell", "Hypodermis",
               "Dermal mesenchyme",
               "Periderm", "Suprabasal Cell")

# filter df
hypoth_filt_df <- dea_ligand_res %>% 
  dplyr::filter(genotype == "hypoth") %>%
  dplyr::filter(cell_type %in% celltypes)

hypo_mat <- table(hypoth_filt_df$pathway, hypoth_filt_df$cell_type)

pdf("plots/fig6b_hypo_lig-wt-up_heatmap.pdf", width = 5.5, height = 4)
pheatmap(hypo_mat, border_color = NULL, cluster_rows = F, 
         color = colorRampPalette(colors = hm_colors)(50))
dev.off()


# Figure 6c ---------------------------------------------------------------

##  plot DE expression of TH effectors ##

# subset cds by cells and genes
sub_cds <- all_cds[,colData(all_cds)$genotype %in% c("wt", "hypoth") & 
                    colData(all_cds)$cell_type_broad == "Basal Cell"]

ncol(sub_cds) # 14k cells
colData(sub_cds)$genotype <- factor(colData(sub_cds)$genotype, levels = c("wt", "nkt", "hypoth"))
my_genes <- get.gene.ids(sub_cds, c("pdgfaa", "pdgfba"))
cds_subset <- sub_cds[my_genes, ]

# plot
plot_percent_cells_positive(cds_subset, 
                                 group_cells_by = "genotype", nrow = 1, ncol = 3) +
  scale_fill_manual(values = geno.colors) +
  theme(legend.position = "none", axis.title.x = element_blank())

ggsave("plots/fig6c_wt-hypo_basal-cell_pdgfaa-ab_perc-positive.pdf", 
       width = 3, height = 2)


# Figure 6d ---------------------------------------------------------------

# filter cell types and clean a duplicate gene
subset_cds <- wt_cds[,partitions(wt_cds) %in% c(1, 4)]
clean_cds <- subset_cds[rowData(subset_cds)$id != "ENSDARG00000110069",]

# load custom plotting function 
plot.pdgf.genes = function(cds, gene1, gene2, gene3, gene4, thresh = 1) {
  gene.id1 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene1, "id"])
  gene.id2 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene2, "id"]) 
  gene.id3 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene3, "id"]) 
  gene.id4 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene4, "id"]) 
  
  colData(cds)$umap_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$umap_2 <- reducedDims(cds)[["UMAP"]][,2]
  
  coldat = colData(cds) %>% 
    as.data.frame()
  
  coldat = coldat %>% 
    mutate(tmp = case_when(assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh ~ "both_lig",
                           assay(cds)[gene.id1,] >= thresh ~ gene1,
                           assay(cds)[gene.id2,] >= thresh ~ gene2,
                           assay(cds)[gene.id3,] >= thresh & assay(cds)[gene.id4,] >= thresh ~ "both_rec",
                           assay(cds)[gene.id3,] >= thresh ~ gene3,
                           assay(cds)[gene.id4,] >= thresh ~ gene4,
                           TRUE ~ "none"))
  
  print(table(coldat$tmp))
  
  coldat$tmp = factor(coldat$tmp, levels = c(gene1, gene2, gene3, gene4, "both_rec", "both_lig", "none"))
  
  plot = ggplot() +
    geom_point(data = coldat %>% 
                 filter(tmp == "none"),
               aes(x = umap_1, y = umap_2), color = "grey80", 
               size = 0.5, stroke = 0,
               alpha = 0.5) +
    geom_point(data = coldat %>% 
                 filter(tmp != "none"),
               aes(x = umap_1, y = umap_2, color = tmp),
               size = 0.5, stroke = 0, alpha = 0.9) +
    scale_color_manual(values = c("#C2A5CF", "#C2A5CF", '#FEDA8B', '#FEDA8B', "#C2A5CF", '#FEDA8B')) +
    #scale_alpha_manual(values = c(1, .5)) + guides(alpha = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    monocle3:::monocle_theme_opts() +
    theme_void()
  
  return(plot)
}

# plot
plot.pdgf.genes(clean_cds, "pdgfaa", "pdgfab", "pdgfra", "pdgfrb", 
                thresh = 1) +
  theme(legend.position = "none")

ggsave("plots/fig6d_wt_pdgf_lig-rec_2color_umap.png", width = 1.5, height = 2,
       dpi = 750, bg = "transparent")


