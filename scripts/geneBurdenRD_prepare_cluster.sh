#!/bin/sh
#$ -o cluster/geneBurdenRD_prepare_$JOB_ID.stdout
#$ -e cluster/geneBurdenRD_prepare_$JOB_ID.stderr
#$ -cwd
#$ -pe smp 1       # request n CPU core
#$ -l h_vmem=1G    # request x Gb RAM / core
#$ -l h_rt=1:0:0  # request x hours runtime

module load R/4.1.0

# Specify 2 mandatory R argument in quotes:
# arg[1]. path for Exomiser/Genomiser candidate variant file (all cases and controls):
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/exomiser_master_file_passvars.tsv"
# arg[2]. path for final output file which will be the input data file for matrix script:
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/exomiserPassWide.tsv"

Rscript scripts/geneBurdenRD_prepare.R "data/exomiser_master_file_passvars.tsv" "data/exomiserPassWide.tsv"

# To run this script on a cluster, from 'geneBurdenRD' folder:
# qsub scripts/geneBurdenRD_prepare_cluster.sh
