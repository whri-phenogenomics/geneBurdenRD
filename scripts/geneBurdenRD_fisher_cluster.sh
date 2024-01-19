#!/bin/sh
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -pe smp 1       # request n CPU core
#$ -l h_vmem=1G    # request x Gb RAM / core
#$ -l h_rt=1:0:0  # request x hours runtime

# Specify task from 1 to total number of different case-control analyses (as per in 'analysisLabelList' file)
#$ -t 1

scriptname="geneBurdenRD_fisher"
exec >cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stdout 2>cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stderr

module load R/4.1.0

# Specify 2 mandatory R arguments:
# args[1]. job index: e.g. $SGE_TASK_ID
# arg[2]. path for file containing list of case-control analysis labels in quotes: e.g. "data/analysisLabelList.tsv"

Rscript scripts/geneBurdenRD_fisher.R $SGE_TASK_ID "data/analysisLabelList.tsv"

# To run this script on a cluster, from 'geneBurdenRD' folder:
# qsub scripts/geneBurdenRD_fisher_cluster.sh
