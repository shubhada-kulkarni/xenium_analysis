#!/bin/bash
#SBATCH --job-name="Xenium_processing"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"

# arguments: path  to Xenium folder (unzipped)

cd /prj/XeniumProbeDesign/analysis_scripts/

module load R/4.4.1_deb12

Rscript --vanilla process_Xenium_data.R $1
