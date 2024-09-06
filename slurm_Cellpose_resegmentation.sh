#!/bin/bash
#SBATCH --job-name="Cellpose_segmentation"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
# #SBATCH --gres=gpu:pascal:1
# #SBATCH --partition=gpu
#SBATCH --mem=100G
# #SBATCH --nodelist=gpu-g2-1

# arguments
# 1. folder with Xenium output for which the segmentation needs to be done

in_folder=$1    # xenium folder

# first go to the folder where input data is stored
cd $in_folder

date
## training
# time  python -m cellpose --dir $in_folder --train --chan 0 --chan2 0 --diameter 10 --do_3D --save_tif --verbose --img_filter _morphology.ome --save_tif --use_gpu --mask_filter "_seg.npy"
echo "Starting cell segmentation"
time python -m cellpose --dir $in_folder --pretrained_model cyto --chan 2 --chan2 0 --diameter 5 --do_3D --save_tif --verbose --img_filter _0001.ome

# after this command is successful, also run the map transcript function to create cell-feature matrix
echo "Transcript mapping"
date
 
#python /prj/XeniumProbeDesign/analysis_scripts/map_transcripts_Cellpose.py -cellpose cellpose_level_3_morphology.ome_seg.npy -transcript transcripts.parquet -out feature_cell_matrix -pix_size 1.7
date
echo "Done!"
