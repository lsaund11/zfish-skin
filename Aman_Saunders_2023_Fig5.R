## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to generate plots for Figure 5 from 
## processed data files available in the GEO repository GSE224695.

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(forcats)
})

source("Aman_Saunders_2023_utils.R")

# set genotype colors
geno.colors <- c("wt" = "#0077BB",
                 "nkt" = "#DDAA33",
                 "hypoth" = "#BB5566",
                 "bon" = "#009988")

# Load data ---------------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")

# Figure 5a ---------------------------------------------------------------

## plot eda and edar in epidermal and dermal cells ##

sub_cds <- all_cds[,colData(all_cds)$genotype %in% c("wt", "bon")] # note: bon epidermal/dermal cells are effectively wildtype
sub_cds <- sub_cds[,partitions(sub_cds) %in% c(5, 3, 2)]

ncol(sub_cds) # 32k cells

# load custom plotting function
plot_fig4_genes = function(cds, gene1, gene2, gene3, thresh = 1, cell_size = 0.5) {
  gene.id1 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene1, "id"])
  gene.id2 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene2, "id"])
  gene.id3 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene3, "id"])
  
  colData(cds)$umap_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$umap_2 <- reducedDims(cds)[["UMAP"]][,2]
  
  coldat = colData(cds) %>% 
    as.data.frame()
  
  coldat = coldat %>% 
    mutate(expr = case_when(assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh ~ "gene1+gene2",
                            assay(cds)[gene.id1,] >= thresh ~ gene1,
                            assay(cds)[gene.id2,] >= thresh ~ gene2,
                            assay(cds)[gene.id3,] >= thresh ~ gene3,
                            TRUE ~ "none"))
  
  coldat$expr = factor(coldat$expr, levels = c(gene1, "gene1+gene2", gene2, gene3, "none"))
  
  plot = ggplot() +
    geom_point(data = coldat %>%
                 filter(expr == "none"),
               aes(x = umap_1, y = umap_2), color = "grey90", 
               size = cell_size, stroke = 0,
               alpha = 1) +
    geom_point(data = coldat %>% 
                 dplyr::arrange(forcats::fct_rev(expr)) %>% 
                 dplyr::filter(expr != "none"),
               aes(x = umap_1, y = umap_2, color = expr),
               size = cell_size, stroke = 0, alpha = 0.7) +
    scale_color_manual(values = c("#ef233c", "#ef233c", "#4895ef", "#52b788")) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    monocle3:::monocle_theme_opts() +
    theme_void()
  
  return(plot)
}

# plot expression of 3 genes
plot_fig4_genes(sub_cds, "shha", "edar", "eda", thresh = 1.1, cell_size = 0.4) +
  theme(legend.position = "none")

ggsave("plots/fig5a_wt_epi_derm_eda-edar-shha_expr_umap.png",
       dpi = 500, width = 1.5, height = 2.2, bg = "transparent")


# Figure 5b ---------------------------------------------------------------

## FGF ligand signature in wildtype vs. eda basal cells ##

# 1. make cell subset
sub_cds <- all_cds[,colData(all_cds)$genotype %in% c("wt", "nkt")]
sub_cds <- sub_cds[,partitions(sub_cds) %in% c(5, 3, 2)]

# 1. identify fgf ligands that are preferentially expressed in basal cells 
fgf_ligs_all <- rowData(all_cds) %>% 
  as.data.frame() %>% 
  filter(grepl("fgf[[:digit:]]", gene_short_name)) %>%
  pull(id)

fgf_cds <- all_cds[fgf_ligs_all,]

# 2. calculate specificity scores for all fgf ligands in each cell type
marker_test_res <- top_markers(fgf_cds, 
                              genes_to_test_per_group = 32, # this is all of them
                              group_cells_by = "cell_type_broad")

bc_fgfs <- marker_test_res %>% 
  filter(cell_group == "Basal Cell") %>% 
  filter(specificity > 0.1) %>% 
  select(gene_short_name) %>% 
  mutate(group = "1")

# duplicate groups to appease aggregation function
bc_fgfs <- rbind(bc_fgfs, bc_fgfs %>% 
                   mutate(group = "2"))

# 3. make signature
fgf_aggr_expr <- aggregate_gene_expression(sub_cds, 
                                          gene_group_df = bc_fgfs, 
                                          scale_agg_values = FALSE)

colData(sub_cds)$fgf_signature <- fgf_aggr_expr[1,]

# extract metadata and plot
coldata_df <- colData(sub_cds) %>% 
  as.data.frame()  %>% 
  filter(genotype %in% c("wt", "nkt"))


ggplot(coldata_df, aes(umap1, umap2, color = fgf_signature)) +
  geom_point(size = 0.7, stroke = 0, alpha = 1, color = "grey90") +
  geom_point(data = subset(coldata_df, fgf_signature > 0.01), 
             aes(x = umap1, y = umap2, color = fgf_signature), size = 0.7, stroke = 0) +
  scale_color_viridis_c() +
  xlim(2.5, 10) +
  ylim(-10,-2.5) +
  facet_wrap(~factor(genotype, levels = c("wt", "nkt")), ncol = 1) +
  theme_void() +
  theme(legend.position = "right")

ggsave("plots/fig5b_wt-nkt-hypo_bc_fgf-signature_umap.png", 
       dpi = 600, bg = "transparent", 
       width = 2, height = 3.5)


# Figure 5c ---------------------------------------------------------------

## generate a positive signature barplot ##

# the function requires gene expression as input so we'll 
# make a little fake cds with the signature score as the gene instead
meta_tmp <- colData(sub_cds) %>% 
  as.data.frame() 

counts_tmp <- meta_tmp %>% 
  select(fgf_signature)
counts_tmp <- t(as.matrix(counts_tmp))
counts_tmp[,1:5]

genes_tmp <- data.frame(id = rownames(counts_tmp), 
                        gene_short_name = rownames(counts_tmp))
rownames(genes_tmp) <- genes_tmp$id

cds_tmp <- new_cell_data_set(expression_data = counts_tmp,
                             gene_metadata = genes_tmp,
                             cell_metadata = meta_tmp)

colData(cds_tmp)$genotype <- factor(colData(sub_cds)$genotype, levels = c("wt", "nkt"))

# plot
plot_percent_cells_positive(cds_tmp, 
                                 group_cells_by = "genotype", normalize = F) +
  scale_fill_manual(values = geno.colors) +
  theme(legend.position = "none", axis.title.x = element_blank())

ggsave("plots/wt-nkt_basal-cell_fgf-signature_perc-positive.pdf", 
       width = 1.5, height = 2)

