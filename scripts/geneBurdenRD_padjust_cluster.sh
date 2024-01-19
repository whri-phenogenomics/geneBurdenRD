#!/bin/sh
#$ -o cluster/geneBurdenRD_padjust_$JOB_ID.stdout
#$ -e cluster/geneBurdenRD_padjust_$JOB_ID.stderr
#$ -cwd
#$ -pe smp 1       # request n CPU core
#$ -l h_vmem=1G    # request x Gb RAM / core
#$ -l h_rt=1:0:0  # request x hours runtime

module load R/4.1.0

# Specify 1 mandatory R argument in quotes:
# args[1]. project name (e.g. same name of project folder): "geneBurdenRD"

Rscript scripts/geneBurdenRD_padjust.R "geneBurdenRD"

# To run this script on a cluster, from 'geneBurdenRD' folder:
# qsub scripts/geneBurdenRD_padjust_cluster.sh
