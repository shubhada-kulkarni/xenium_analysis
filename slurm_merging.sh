#!/bin/bash
#SBATCH --job-name="Xenium_merging"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=shubhada.kulkarni@uni-heidelberg.de
#SBATCH --error="%x.err.%j"
#SBATCH --output="%x.out.%j"

# arguments:
# 1. region1 xenium processed object
# 2. region2 xenium processed object
# 3. slide name for output file

cd /prj/XeniumProbeDesign/analysis_scripts/

module load R/4.4.1_deb12

Rscript --vanilla merge_regions.R $1 $2 $3
