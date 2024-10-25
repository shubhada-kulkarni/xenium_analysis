# R script to process Xenium data stored in given path

library(Seurat, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")

args = commandArgs(trailingOnly = TRUE)
print(args)
path <- args[1]
id <- args[2]
print(id)

# storing RDS object of loaded Xenium object
# id <- basename(path)
out_path <- "/prj/XeniumProbeDesign/heart_human_29072024/xenium_objects_R/"
outRDS <- paste(out_path, id, "_full.rds", sep="")
print(outRDS)

# read in the Xenium data
print("Loading in the Xenium data")
xenium.obj <- LoadXenium(path, fov = "fov")
saveRDS(object = xenium.obj, file = outRDS)

# subsetting the Xenium data to remove low counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 10) # this filtering needs to be done or there is an error

# Processing
print("Processing the Xenium data")
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj)

# saving processed object
print("Saving the processed Xenium object")
outRDS_processed <- paste(out_path, id, "_processed.rds", sep="")
saveRDS(xenium.obj, file = outRDS_processed)

# colouring cells by clusters
ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)

print("Done!")
