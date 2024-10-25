#!/bin/bash
#SBATCH --job-name="Xenium_import"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"

# arguments
# 1. folder with Xenium output for which the segmentation needs to be done

in_folder=$1    # xenium folder
gson_nuclei=$2  # GeoJSON polygons for nuclei
gson_cells=$3   # GeoJSON polygons for cell boundaries

## first go to the folder where input data is stored
cd $in_folder

xeniumranger import-segmentation --xenium-bundle=$in_folder --id=qupath_seg --nuclei=$gson_nuclei --cells=$gson_cells

##########################
#
# cd /prj/XeniumProbeDesign/analysis_scripts

# xeniumranger import-segmentation --xenium-bundle=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/ --id=qupath --nuclei=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/morphology_focus_0000.geojson --cells=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/morphology_focus_0003.geojson
