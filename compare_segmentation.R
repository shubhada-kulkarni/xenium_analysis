# r script to compare results from original segmentation of Xenium and re-segmentation using Baysor
library(Seurat)
library(R.utils)

original <- LoadXenium("/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/", fov = "fov")
baysor <- LoadXenium("/prj/XeniumProbeDesign/analysis_scripts/baysor-demo/outs/", fov = "fov")

# processing
baysor <- subset(baysor, subset = nCount_Xenium > 50)
baysor <- SCTransform(baysor, assay = "Xenium")
baysor <- RunPCA(baysor, npcs = 30, features = rownames(baysor))
baysor <- RunUMAP(baysor, dims = 1:30)
baysor <- FindNeighbors(baysor, reduction = "pca", dims = 1:30)
baysor <- FindClusters(baysor, resolution = 0.25)

transcripts_original <- read.csv(gzfile("/prj/XeniumProbeDesign/heart_human_29072024/resegmented/output-XETG00046__0018072__Region_1__20240725__112631_resegmented/outs_0018072_Region1/transcripts.csv.gz"))
transcripts_baysor <- read.csv(gzfile("/prj/XeniumProbeDesign/analysis_scripts/baysor-demo/outs/transcripts.csv.gz"))

# check what exactly is cluster 5 by its markers
markers_5 <- FindMarkers(baysor, ident.1 = "5")

