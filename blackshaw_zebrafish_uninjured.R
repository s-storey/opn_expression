library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(data.table)
library(magrittr)
library(Matrix)


# dataset downloaded from https://github.com/jiewwwang/Single-cell-retinal-regeneration
# sample information 
  # Zebrafish_gene_features.tsv 
  # Zebrafish_LD_cell_features.tsv 
  
# Count matrix --> http://bioinfo.wilmer.jhu.edu/jiewang/scRNAseq/Count_matrix/
  # Zebrafish_LD_count_matrix.mtx

fish <- ReadMtx(
  mtx = 'Zebrafish_LD_count_matrix.mtx',
  cells = 'Zebrafish_LD_cell_features.tsv',
  features = 'Zebrafish_gene_features.tsv',
  cell.column = 2,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 1,
  skip.feature = 1,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)

fish_obj = CreateSeuratObject(counts = fish)

metadata <- read.csv('meta.csv') ## adapted from Zebrafish_LD_cell_features.tsv
Barcode <- metadata[1]
Sample <- metadata[2]
Cell.type<- metadata[3]
tSNE_1 <- metadata[4]
tSNE_2 <- metadata[5]
nGene <- metadata[6]
nUMI <- metadata[7]
Percentage.of.mitochondrial.genes <- metadata[8]
Percentage.of.ribosomal.protein.genes <- metadata[9]


fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(Barcode = Barcode, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(Sample = Sample, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(Cell.type = Cell.type, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(tSNE_1 = tSNE.1, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(tSNE_2 = tSNE.2, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(nGene = nGene, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(nUMI = nUMI, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(Percentage.of.mitochondrial.genes = Percentage.of.mitochondrial.genes, row.names = rownames(fish_obj@meta.data)))
fish_obj <- AddMetaData(object = fish_obj, metadata = data.frame(Percentage.of.mitochondrial.genes = Percentage.of.mitochondrial.genes, row.names = rownames(fish_obj@meta.data)))


Idents(object = fish_obj) <- "Sample"

control_fish <- subset(x = fish_obj, idents = c("Adult R1","Adult R2", "Adult R3", "Adult R4", "Adult R5" ))

qc_fish <- subset(control_fish, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & Percentage.of.mitochondrial.genes < 15)

subset.list <- SplitObject(qc_fish, split.by = 'Sample')
for (i in 1:length(subset.list)) {
  subset.list[[i]] <- NormalizeData(subset.list[[i]], verbose = FALSE) 
  subset.list[[i]] <- FindVariableFeatures(subset.list[[i]], selection.method = "vst",
                                           nfeatures = 2000, verbose = FALSE)
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = subset.list)

integration.anchors <- FindIntegrationAnchors(object.list = subset.list, anchor.features = features)
# this command creates an 'integrated' data assay
integrated.combined <- IntegrateData(anchorset = integration.anchors)

saveRDS(integrated.combined, 'integrated.rds')

int_fish <- integrated.combined


int_fish

DefaultAssay(int_fish) <- 'integrated'


all.genes <- rownames(int_fish)
int_fish <- ScaleData(int_fish, features = all.genes)
int_fish <- RunPCA(int_fish)
int_fish <- FindNeighbors(int_fish)
int_fish <- FindClusters(object = int_fish)
int_fish <- RunTSNE(int_fish)
int_fish <- RunUMAP(int_fish, reduction = "pca", dims = 1:19)
DefaultAssay(int_fish) <- 'RNA'

# Set the active identity to "Cell.type"
int_fish <- SetIdent(int_fish, value = int_fish@meta.data$Cell.type)
DimPlot(object = int_fish, reduction = "umap", label = TRUE)

saveRDS(int_fish, 'control_fish_integrated.rds')