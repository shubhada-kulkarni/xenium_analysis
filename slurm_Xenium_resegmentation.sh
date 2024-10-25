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

# /prj/XeniumProbeDesign/analyser_installation/xeniumranger-xenium2.0/xeniumranger resegment --xenium-bundle $in_folder --id $out_id
##########################
# Do not use boundary stain in analysis, but keep default interior stain and DAPI
/prj/XeniumProbeDesign/analyser_installation/xeniumranger-xenium2.0/xeniumranger resegment --xenium-bundle $in_folder --id $out_id --boundary-stain=disable --interior-stain=18S --localcores=64 --localmem=256

##########################

## Do not use either boundary or interior stain and only use DAPI for nuclear expansion
# (this option worked the best till now)
#/prj/XeniumProbeDesign/analyser_installation/xeniumranger-xenium2.0/xeniumranger resegment --xenium-bundle $in_folder --id $out_id --boundary-stain=disable --interior-stain=disable --localcores=128 --localmem=256

## trying above option with larger expansion distances
#/prj/XeniumProbeDesign/analyser_installation/xeniumranger-xenium2.0/xeniumranger resegment --xenium-bundle $in_folder --id $out_id --boundary-stain=disable --interior-stain=disable --localcores=128 --localmem=256 --expansion-distance=10
##########################
