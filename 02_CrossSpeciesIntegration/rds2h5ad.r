library(rhdf5)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
prefix <- args[2]


seu.obj <- readRDS(rds_file)

#if (startsWith(as.character(seu.obj@version), "5")) {
    #seu.obj <-  JoinLayers(seu.obj)
#}

rna_assay <- GetAssay(seu.obj, "RNA")


counts <- GetAssayData(rna_assay, "counts")


mtx_file <- paste(prefix, "mtx", sep=".")
row_file <- paste(prefix, "rows", sep=".")
col_file <- paste(prefix, "cols", sep=".")
metadata_file <- paste(prefix, "metadata.csv", sep=".")


# write MM file for Seurat
Matrix::writeMM(counts, mtx_file)
# write rowname and colnames
write.table(rownames(counts), row_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(colnames(counts), col_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
write.csv(seu.obj[[]], metadata_file)
