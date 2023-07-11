library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(data.table)

# This function will generate a dataframe including annotated clusters, number of cells
# in each cluster, and the number of cells expressing each gene from a list in the cluster.


# Instructions: import suerat object from .rds file
seurat_object <- readRDS('your_file_here.rds')


# Create list of genes ('master' in example below)


master <- c('rho', 'rhol', 'opn1sw1', 'opn1sw2', 'opn1mw1', 'opn1mw2', 'opn1mw3', 'opn1mw4', 'opn1lw1', 'opn1lw2', 'opn4.1', 'opn4a', 'opn4b','opn5', 'opn6a', 'opn6b',  'opn7c', 'opn7d', 'opn8a', 'opn8b', 'opn8c', 'opn9')


############################## IMPORTANT #######################################
#### This funciton will create a dataframe of annotated clusters            ####
#### Annotation in this example seurat object is 'CellType'.                ####
#### Please change 'CellType' in the below function to match your           ####
#### annotation.                                                            ####
################################################################################
fish_opn_expression <- function(seurat_object, gene_list, min_expression) {
  # Get cell counts by cluster
  counts <- as.data.frame(table(seurat_object@meta.data$CellType)) # <-- change 'CellType' to appropriate metadata column for annotation
  
  # Initialize a matrix to store the gene expression counts
  gene_counts <- matrix(NA, nrow = nrow(counts), ncol = length(gene_list), 
                        dimnames = list(rownames(counts), gene_list))
  
  # Calculate gene expression counts
  for (i in 1:length(gene_list)) {
    gene <- gene_list[i]
    tryCatch({
      gene_counts[, i] <- apply(counts, 1, function(x) {
        cluster <- x[1]
        num_expressing_cells <- sum(seurat_object@assays$RNA@counts[gene, seurat_object@meta.data$CellType == cluster] >= min_expression)# <-- change 'CellType' to appropriate metadata column for annotation
        return(num_expressing_cells)
      })
    }, error = function(e) {
      message(paste("Error: Failed to process gene", gene, ". Skipping."))
    })
  }
  
  # Combine gene expression counts with cluster counts
  counts <- cbind(counts, gene_counts)
  
  # Rename columns
  colnames(counts)[-c(1, 2)] <- paste0("Cells_expressing_", gene_list)
  
  return(counts)
}

# Generate dataframe by running function with appropriate object, gene list and
# minimum expression level. 

gene_counts_df <- fish_opn_expression(blackshaw_fish, master, 1)

View(gene_counts_df)

