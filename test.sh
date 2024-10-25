#!/bin/bash
#
#SBATCH --job-name="baysor"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
#SBATCH --mem=100G

echo "Running baysor import stains"
date

cd /prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/
mkdir test_imageJ_0.5

# ~/.julia/bin/baysor run formatted_transcripts.csv C2-QuPath_seg_boundary.tif -x x_location -y y_location -z z_location -g feature_name -m 10 --prior-segmentation-confidence 0.8 --save-polygons GeoJSON -o test_imageJ_0.8
~/.julia/bin/baysor run formatted_transcripts.csv C2-QuPath_seg_boundary.tif -x x_location -y y_location -z z_location -g feature_name -m 10 --prior-segmentation-confidence 0.5 --save-polygons GeoJSON -o test_imageJ_0.5


date
echo "Done"

## removed commands
# xeniumranger import-segmentation --id=cellpose-seg --xenium-bundle=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/ --cells=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/cellpose_level_3_morphology.ome_seg.npy --nuclei=/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/morphology_focus/morphology_focus_0000.ome.tif --localcores=32 --localmem=128
#
#
#
# ~/.julia/bin/baysor run formatted_transcripts.csv RGB_seg.tif -x x_location -y y_location -z z_location -g feature_name -m 10 --prior-segmentation-confidence 0.8 --save-polygons GeoJSON -o test
