# R script to read in two regions from one slide and merge the data

options(future.globals.maxSize = 8000 * 1024^2) # for memory reasons

library(Seurat, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")
library(Rfast2, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")
library(dplyr, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")

args = commandArgs(trailingOnly=TRUE)
region1_file <- args[1] #"/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/output-XETG00046__0027119__Region_1__20240621__120943_processed.rds")
region2_file <- args[2] #"/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/output-XETG00046__0027119__Region_2__20240621__120943_processed.rds"
slide_name <- args[3]

print("Reading individual region objects")
region1 <- readRDS(region1_file)
region2 <- readRDS(region2_file)

print("Merging")
xenium.obj <- merge(region1, y = region2, add.cell.ids = c("region1", "region2"), project = "slide_name")

out_rds <- paste("/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/", slide_name, "_merged.rds", sep = "")
saveRDS(xenium.obj, out_rds)

print("Process the merged data")

xenium.obj@meta.data$region <- unlist(lapply(rownames(xenium.obj@meta.data), function(x) strsplit(x, "_")[[1]][1]))
features <- Features(xenium.obj)
print("Normalise and scale the merged data")
#xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- ScaleData(xenium.obj)
# variable features set to all genes,  for PCA calculations
VariableFeatures(xenium.obj) <- features
print("PCA and UMAP calculations")
xenium.obj <- RunPCA(xenium.obj)
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
# Plot to check if the regions perfectly overlap.  If they do, there is no need of batch effect correction
png(paste("/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/slide_", slide_name, "_dimplot_regions_overlap.png", sep = ""), width = 800, height = 800)
DimPlot(xenium.obj, group.by = "region")  
dev.off()
print("Finding neighbors and clustering analysis")
xenium.obj <- FindNeighbors(xenium.obj, dims = 1:30) 
xenium.obj <- FindClusters(xenium.obj)

print("Saving processed object")
out_rds_processed <- paste("/prj/XeniumProbeDesign/kidney_Nphs2-mice_Xenium_Martin/xenium_objects_R/", slide_name, "_merged_processed.rds", sep = "")
saveRDS(xenium.obj, out_rds_processed)

print("DONE!")
