#!/bin/bash
#
#SBATCH --job-name="baysor"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"
#SBATCH --mem=100G

cd /prj/XeniumProbeDesign/analysis_scripts/

file_transcript=$1
out_formatted="formatted_"$file_transcript

echo "formatting the CSV file"
# first format the input transcript csv file (unzipped)
python filter_transcripts_Baysor.py -transcript $file_transcript -out $out_formatted

echo "running baysor"
# next supply this formatted output to baysor run
~/.julia/bin/baysor run -x x_location -y y_location -z z_location -g feature_name -m 10 -p --prior-segmentation-confidence 0.1 $out_formatted :cell_id_int --save-polygons GeoJSON -o /prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/baysor/