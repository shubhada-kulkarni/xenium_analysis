#!/bin/bash
#SBATCH --job-name="Xenium_segmentation"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"

# arguments
# 1. folder with Xenium output for which the segmentation needs to be done

in_folder=$1    # xenium folder
out_id=$2       # output re-segmented folder (without full path)

# first go to the folder where input data is stored
path=`dirname $in_folder`

cd $path

/prj/XeniumProbeDesign/analyser_installation/xeniumranger-xenium2.0/xeniumranger resegment --xenium-bundle $in_folder --id $out_id
