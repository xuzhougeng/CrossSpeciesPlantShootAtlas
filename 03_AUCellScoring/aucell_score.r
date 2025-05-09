#BiocManager::install("AUCell")
library(Seurat)
library(AUCell)
library(GSEABase)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript aucell_score.r <seurat_rds> <geneset_csv>")
}

seurat_file <- args[1]
geneset_file <- args[2]

# Read Seurat object
seu.obj <- readRDS(seurat_file)

# Read geneset from CSV
geneset_df <- read.csv(geneset_file, header = TRUE)
geneset_list <- split(geneset_df$marker_gene, geneset_df$cell_type)

# Process each cell type
for(cell_type in names(geneset_list)) {
  genes <- unique(geneset_list[[cell_type]])
  geneSets <- GeneSet(genes, setName=cell_type)
  
  # Calculate scores
  assay <- GetAssay(seu.obj, "RNA")
  exprMatrix <- GetAssayData(assay, "counts")
  cells_AUC <- AUCell_run(exprMatrix, geneSets)
  seu.obj[[paste0("score_", cell_type)]] <- getAUC(cells_AUC)[cell_type,]
  
  # Save plot
  png(paste0("score_", cell_type, ".png"), width=800, height=600)
  print(FeaturePlot(seu.obj, paste0("score_", cell_type)))
  dev.off()
}