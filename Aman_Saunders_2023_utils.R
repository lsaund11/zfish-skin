# Utility functions for Saunders, Srivatsan, et al., 2023

# extract coordinates for a PAGA graph from a cell data set 
get_paga_graph <- function(cds, reduction_method = "UMAP") {
  
  cluster_result <- cds@clusters[[reduction_method]]$cluster_result
  
  cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                     cluster_result$optim_res,
                                                     qval_thresh = 0.05, 
                                                     FALSE)
  
  cluster_g <- cluster_graph_res$cluster_g
  cluster_g <- igraph::set_vertex_attr(cluster_g, "name", value = stringr::str_replace(igraph::V(cluster_g)$name, "cell_membership", ""))
  cluster_g
}


# add 2d or 3d umap coordinates to your coldata
append_umap_coordinates = function(cds, umap_3D = F){
  
  if (!umap_3D){
    colData(cds)$umap1 = reducedDim(x = cds,
                                    type = "UMAP")[,1]
    colData(cds)$umap2 = reducedDim(x = cds,
                                    type = "UMAP")[,2]
  }
  if (umap_3D){
    colData(cds)$umap3d_1 = reducedDim(x = cds,
                                       type = "UMAP")[,1]
    colData(cds)$umap3d_2 = reducedDim(x = cds,
                                       type = "UMAP")[,2]
    colData(cds)$umap3d_3 = reducedDim(x = cds,
                                       type = "UMAP")[,3]
  }
  return(cds)
}

# calculate proliferation index based on cell cycle gene expression
estimate_cc_scores = function(cds, gene_list1, 
                              gene_list2){
  
  cds_g1s = cds[fData(cds)$id %in% gene_list1,]
  aggregate_g1s_expression = assay(cds_g1s)
  aggregate_g1s_expression = t(t(aggregate_g1s_expression) / colData(cds_g1s)$Size_Factor)
  aggregate_g1s_expression = Matrix::colSums(aggregate_g1s_expression)
  
  cds_g2m = cds[fData(cds)$id %in% gene_list2,]
  aggregate_g2m_expression = assay(cds_g2m)
  aggregate_g2m_expression = t(t(aggregate_g2m_expression) / pData(cds_g2m)$Size_Factor)
  aggregate_g2m_expression = Matrix::colSums(aggregate_g2m_expression)
  
  pData(cds)$g1s_score = log(aggregate_g1s_expression+1)
  pData(cds)$g2m_score = log(aggregate_g2m_expression+1)
  
  pData(cds)$proliferation_index = log(aggregate_g1s_expression + aggregate_g2m_expression + 1)
  
  return(cds)
  
}

# an updated function for plotting two genes in one plot
plot_two_genes = function(cds, gene1, gene2, thresh = 1, cell_size = 0.5) {
  gene.id1 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene1, "id"])
  gene.id2 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene2, "id"])
  
  colData(cds)$umap_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$umap_2 <- reducedDims(cds)[["UMAP"]][,2]
  
  coldat = colData(cds) %>% 
    as.data.frame()
  
  coldat = coldat %>% 
    mutate(expr = case_when(assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh ~ "both",
                            assay(cds)[gene.id1,] >= thresh ~ gene1,
                            assay(cds)[gene.id2,] >= thresh ~ gene2,
                            TRUE ~ "neither"))
  
  coldat$expr = factor(coldat$expr, levels = c(gene1, gene2, "both", "neither"))
  
  plot = ggplot() +
    geom_point(data = coldat %>% 
                 filter(expr == "neither"),
               aes(x = umap_1, y = umap_2), color = "grey90", 
               size = cell_size, stroke = 0,
               alpha = 1) +
    geom_point(data = coldat %>% 
                 filter(expr != "neither"),
               aes(x = umap_1, y = umap_2, color = expr),
               size = cell_size, stroke = 0, alpha = 0.7) +
    scale_color_manual(values = c("#4895ef", "#ef233c", "#9163cb")) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    monocle3:::monocle_theme_opts() +
    theme_void()
  
  return(plot)
}

# return threshold levels for the plot two genes plot
get_two_genes_df = function(cds, gene1, gene2, thresh = 1) {
  gene.id1 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene1, "id"])
  gene.id2 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene2, "id"])
  
  colData(cds)$umap_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$umap_2 <- reducedDims(cds)[["UMAP"]][,2]
  
  coldat = colData(cds) %>% 
    as.data.frame()
  
  coldat = coldat %>% 
    mutate(expr = case_when(assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh ~ "both",
                            assay(cds)[gene.id1,] >= thresh ~ gene1,
                            assay(cds)[gene.id2,] >= thresh ~ gene2,
                            TRUE ~ "neither"))
  
  coldat$expr = factor(coldat$expr, levels = c(gene1, gene2, "both", "neither"))
  
  return(coldat)
}

# return a list of gene ids from short names
get.gene.ids = function (cds, names) {
  return(as.character(rowData(cds)[rowData(cds)$gene_short_name %in% names, "id"]))
}

