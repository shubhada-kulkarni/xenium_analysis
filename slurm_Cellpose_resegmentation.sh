#!/bin/bash
#SBATCH --job-name="Cellpose_segmentation"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
#SBATCH --gres=gpu:ampere:2
#SBATCH --partition=gpu
#SBATCH --mem=100G


# arguments
# 1. folder with Xenium output for which the segmentation needs to be done

in_folder=$1    # xenium folder

# first go to the folder where input data is stored
path=`dirname $in_folder`

cd $path

time  python -m cellpose --dir $in_folder --train --chan 0 --chan2 0 --diameter 10 --do_3D --save_tif --verbose --img_filter morphology.ome --save_tif --use_gpu

