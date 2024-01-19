#!/bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -pe smp 1       # request n CPU core
#$ -l h_vmem=1G    # request x Gb RAM / core
#$ -l h_rt=1:0:0  # request x hours runtime

# Specify task from 1 to total number of case-control analyses (as per in 'analysisLabelList.tsv' file)
#$ -t 1

scriptname="geneBurdenRD_visualisation"
exec >cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stdout 2>cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stderr

module load R/4.1.0

# Specify 3 mandatory R arguments:
# args[1]. job index: e.g. $SGE_TASK_ID
# args[2]. path for Exomiser masterfile file in quotes: e.g. "data/exomiser_master_file_passvars.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/exomiser_master_file_passvars.tsv"
# args[3]. path for file containing list of case-control analysis labels in quotes: e.g. "data/analysisLabelList.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/analysisLabelList.tsv"

# Two optional arguments
# args[4]. choice of FDR threshold: "0.01" - default is "0.05"
# args[5]. choice of "plot" or "table_only" - default is "plot"

Rscript scripts/geneBurdenRD_visualisation.R $SGE_TASK_ID "data/exomiser_master_file_passvars.tsv" "data/analysisLabelList.tsv" "0.05" "plot"

# To run this script on a cluster, from 'geneBurdenRD' folder:
# qsub scripts/geneBurdenRD_visualisation_cluster.sh
