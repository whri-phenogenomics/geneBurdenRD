#!/bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -pe smp 1       # request n CPU core
#$ -l h_vmem=1G    # request x Gb RAM / core
#$ -l h_rt=1:0:0  # request x hours runtime

# Specify task from 1 to total number of case-control analyses (as per in 'analysisLabelList.tsv' file)
#$ -t 1

scriptname="geneBurdenRD_matrix"
exec >cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stdout 2>cluster/${scriptname}_${JOB_ID}_${SGE_TASK_ID}.stderr

module load R/4.1.0

echo ${SGE_TASK_ID}

# Specify 3 mandatory R arguments:
# args[1]. job index: e.g. $SGE_TASK_ID
# args[2]. data path for Exomiser/Genomiser file in quotes: e.g. "data/exomiserPassWide.tsv"
# args[3]. data path for file containing list of case-control analysis labels in quotes: e.g. "data/analysisLabelList.tsv"

# One optional argument
# args[4]. whether gene list files are provided: "geneListFileON"

Rscript scripts/geneBurdenRD_matrix.R $SGE_TASK_ID "data/exomiserPassWide.tsv" "data/analysisLabelList.tsv"

# To run this script on a cluster, from 'geneBurdenRD' folder:
# qsub scripts/geneBurdenRD_matrix_cluster.sh
