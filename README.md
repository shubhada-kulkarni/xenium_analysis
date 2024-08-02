Scripts for xenium data analysis


__process_Xenium_data.R__ : R script to read and process Xenium data for one region at a time

__slurm_processing.R__ : Bash script to submit __process_Xenium_data.R__ as a SLURM job

__merge_regions.R__ :  R script to perform merging of two regions on a slide and post-processing the merged data

__slurm_merging.sh__ : Bash script to submit __merge_regions.R__ as a SLURM  job

__postprocessing_Kidney_Xenium.qmd and postprocessing_Kidney_Xenium.html__ : Markdown document and report for cluster annotation analysis for kidney data from Cologne

__figs/__ : Folder with all figures

__postprocessing_humanheart_Xenium.qmd and postprocessing_humanheart_Xenium.html__ : Markdown document and report for cluster annotation analysis for pilot Xenium experiment on human patient heart AV node
