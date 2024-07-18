# R script to process Xenium data stored in given path

library(Seurat, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")

args = commandArgs(trailingOnly=TRUE)
path <- args[1]

# storing RDS object of loaded Xenium object
id <- basename(path)
out_path <- "/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/"
outRDS <- paste(out_path, id, "_full.rds", sep="")

# read in the Xenium data
print("Loading in the Xenium data")
xenium.obj <- LoadXenium(path, fov = "fov")
saveRDS(object = xenium.obj, file = outRDS)

# subsetting the Xenium data to remove low counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# Processing
print("Processing the Xenium data")
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.25)

# saving processed object
print("Saving the processed Xenium object")
outRDS_processed <- paste(out_path, id, "_processed.rds", sep="")
saveRDS(outRDS_processed, file = outRDS_processed)

# colouring cells by clusters
ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)

print("Done!")
