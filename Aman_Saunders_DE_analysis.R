## Aman, Saunders et al. eLife (2023). [https://doi.org/10.7554/eLife.86670.3]
## This script contains code to run the differential
## analysis as described in the paper
## This analysis was run with 250GB of memory
## (5x50GB) on a compute cluster (note: smaller chunks will error out)

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(pheatmap)
  library(devtools)
})

source("Aman_Saunders_2023_utils.R")

# set genotype colors
geno.colors <- c("wt" = "#0077BB",
                "nkt" = "#DDAA33",
                "hypoth" = "#BB5566",
                "bon" = "#009988")


# Load data ---------------------------------------------------------------

all_cds <- readRDS("data/GSM7029635_zskin_all_genotypes_cds.RDS")

# Load gene resource info -------------------------------------------------

# full gene names
gene_names <- read.csv("../bin/zfish_full_gene_names.csv", header = T, stringsAsFactors = F)

# TF and LR lists
tf_lr_df <- read.csv("data/zfish_lig_rec_tf_list.csv", stringsAsFactors = F)

# Wildtype vs. hypoTH DEA (all cell types) ---------------------------------------

# subset and clean data
wt_hypo_sub <- cds[,colData(cds)$genotype %in% c("wt", "hypoth") & !is.na(colData(cds)$cell_type_broad)]
wt_hypo_sub <- estimate_size_factors(wt_hypo_sub)
colData(wt_hypo_sub)$genotype <- factor(colData(wt_hypo_sub)$genotype, levels = c("wt", "hypoth"))
cell.groups <- colData(wt_hypo_sub) %>% 
  as.data.frame() %>% 
  filter(!(cell_type_broad %in% c("Doublet", "Unknown"))) %>% 
  group_by(cell_type_broad) %>% 
  tally() %>% 
  filter(n > 50) %>% 
  pull(cell_type_broad)

# run DE testing
hypoth_cells_list <- list()

for(type in cell.groups){
  message("Finding DEGs for ", type)
  celltype_cds <- wt_hypo_sub[,colData(wt_hypo_sub)$cell_type_broad == type] %>% 
    estimate_size_factors() %>% 
    detect_genes()
  celltype_cds <- celltype_cds[Matrix::rowSums(SingleCellExperiment::counts(celltype_cds) > 0) > 10,]
  fits <- fit_models(celltype_cds, model_formula_str = "~genotype", cores = 1)
  mod_coef <- coefficient_table(fits)
  mod_coef$cell_type <- type
  hypoth_cells_list[[type]] <- mod_coef
}

all.TH.celltype.degs <- data.table::rbindlist(hypoth_cells_list, 
                                             idcol = "cell_type") %>% 
  as.data.frame() %>% 
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::select(-model_summary, -model)

# generate a list of filtered, organized DEGs
all.TH.celltype.degs = all.TH.celltype.degs %>% 
  dplyr::mutate(up_in = case_when(
    estimate < 0 ~ "WT",
    estimate > 0 ~ "hypoTH"))

# save full results
# fwrite(all.TH.celltype.degs, "DEA/wt-hypoth_all-celltype_fit-model_results.csv", sep = ",")

# save DEG counts per celltype
all.TH.celltype.degs %>%
  dplyr::filter(q_value < 0.05) %>% 
  group_by(cell_type) %>% dplyr::summarize(deg_count = n()) %>% arrange(-deg_count) %>% 
  write.csv("DEA/wt-hypo_DEG-q05_celltype_counts.csv", row.names = F)

# merge with gene full names and tf/lig/rec identities
all.TH.celltype.degs %>% 
  filter(term != "(Intercept)") %>% 
  filter(q_value < 0.05) %>%
  arrange(cell_type, up_in, q_value) %>%
  left_join(names, by = "gene_short_name") %>% 
  select(cell_type, up_in, id, gene_short_name, 
         gene_full_name, q_value, estimate) %>%
  left_join(tf_lr_df, by = "id") 

fwrite("DEA/wt_hypo_all-celltypes_DEG_list", sep = ",")

# Tabulate multi-hit DEGs across cell types -------------------------------

# collapse genes that are seen multiple times across celltypes - are there commonalities?
wt.hypo.degs.by.type <- left_join(wt.hypo.multi.df, 
                                  all.TH.celltype.degs %>% 
                                    dplyr::filter(q_value < 0.05) %>%
                                    dplyr::select(id, gene_short_name, cell_type, up_in), 
                                 by = "id") %>% group_by(id) %>% 
  mutate(cell_type = paste(cell_type, collapse=", ")) %>%
  mutate(up_in = ifelse(test = n_distinct(up_in)==1, 
                        yes = up_in[1], no = paste(up_in, collapse = ", "))) %>% 
  distinct(id, .keep_all = TRUE) %>% arrange(-Freq) %>% 
  left_join(names, by = "gene_short_name") %>% 
  select(Freq, id, gene_short_name, gene_full_name, cell_type, up_in)

fwrite(wt.hypo.degs.by.type, 
       file = "DEA/wt-hypo_multi-degs_by_celltype.csv", 
       sep = ",")


# Note: Code for plotting DEG heatmaps are in the Figure 6 notebook


# All genotypes -----------------------------------------------------------

## the code above can be used as a template for all wt vs. other 
## genotype comparisons from the paper
