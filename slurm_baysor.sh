#!/bin/bash
#
#SBATCH --job-name="baysor"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
#SBATCH --mem=100G

xenium_folder=$1

date
cd /prj/XeniumProbeDesign/analysis_scripts # cd to folder where transcript file exists

echo "Unzipping the trancsripts CSV file"
gzip -d $xenium_folder"/transcripts.csv.gz"
echo "Done"
date

echo "formatting the CSV file"
file_transcript=$xenium_folder"/transcripts.csv"
out_formatted=$xenium_folder"/formatted_transcripts.csv"
echo $file_transcript, $out_formatted
# first format the input transcript csv file (unzipped)
python filter_transcripts_Baysor.py -transcript $file_transcript -out $out_formatted
date

echo "running baysor"
# next supply this formatted output to baysor run
baysor_out=$xenium_folder"/baysor"
~/.julia/bin/baysor run -x x_location -y y_location -z z_location -g feature_name -m 10 -p --prior-segmentation-confidence 0.1 $out_formatted :cell_id_int --save-polygons GeoJSON -o $baysor_out
date

echo "Importing segmentations"
cd $xenium_folder
xeniumranger import-segmentation --id=baysor-segmented --xenium-bundle=$xenium_folder --transcript-assignment=$baysor_out"/segmentation.csv" --viz-polygons=$baysor_out"/segmentation_polygons.json" --units=microns --localcores=32 --localmem=128
date
echo "Done"

echo "Re-zipping the transcript CSV file"
gzip $file_transcript