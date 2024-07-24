# R script to read in two regions from one slide and merge the data

options(future.globals.maxSize = 8000 * 1024^2) # for memory reasons

library(Seurat, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")
library(Rfast2, lib.loc = "/beegfs/homes/skulkarni/R/x86_64-pc-linux-gnu-library/4.3")

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
print("DONE!")