#BSUB -q medium
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -o cluster/geneBurdenRD_matrix_jobarray_%J_%I.stdout
#BSUB -e cluster/geneBurdenRD_matrix_jobarray_%J_%I.stderr
#BSUB -cwd /geneBurdenRD
#BSUB -n 5
#BSUB -J geneBurdenRD_matrix_jobarray[1]
#BSUB -R rusage[mem=80000]
#BSUB -M 80000

module load lang/R/4.1.0-foss-2019b

# Specify 3 mandatory R arguments in quotes:
# args[1]. job index: e.g. 'LSB_JOBINDEX'
# args[2]. data path for Exomiser/Genomiser file: e.g. "data/exomiserPassWide.tsv"
# args[3]. data path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"

# One optional argument
# args[4]. whether gene list files are provided: "geneListFileON"

Rscript scripts/geneBurdenRD_matrix_jobarray.R 'LSB_JOBINDEX' "data/exomiserPassWide.tsv" "data/analysisLabelList.tsv" 

# To run this script:

# 1. Uncomment: runMode <- "cluster" in geneBurdenRD_matrix_jobarray.R
# 2. ssh cluster
# 3. bsub < /geneBurdenRD/scripts/geneBurdenRD_matrix_jobarray.sh
# 4. bjobs
