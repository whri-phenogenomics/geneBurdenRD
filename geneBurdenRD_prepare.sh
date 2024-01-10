#BSUB -q short
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -o cluster/geneBurdenRD_prepare_%J.stdout
#BSUB -e cluster/geneBurdenRD_prepare_%J.stderr
#BSUB -cwd /geneBurdenRD
#BSUB -n 2
#BSUB -R rusage[mem=48000]
#BSUB -M 48000

module load lang/R/4.1.0-foss-2019b

# Specify 2 mandatory R argument in quotes:
# arg[1]. path for Exomiser/Genomiser candidate variant file (all cases and controls): 
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/exomiser_master_file_passvars.tsv"
# arg[2]. path for final output file which will be the input data file for matrix script: 
# (path relative to cwd; otherwise, absolute path; in quotes) e.g. "data/exomiserPassWide.tsv"

Rscript /geneBurdenRD/scripts/geneBurdenRD_prepare.R "data/exomiser_master_file_passvars.tsv" "data/exomiserPassWide.tsv"

# To run this script:

# 1. Uncomment: runMode <- "cluster" in geneBurdenRD_prepare.R
# 2. ssh cluster
# 3. /geneBurdenRD/scripts/geneBurdenRD_prepare.sh 
# 4. bjobs 
