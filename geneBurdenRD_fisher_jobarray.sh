#BSUB -q long
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -o cluster/geneBurdenRD_fisher_jobarray_%J_%I.stdout
#BSUB -e cluster/geneBurdenRD_fisher_jobarray_%J_%I.stderr
#BSUB -cwd /re_gecip/enhanced_interpretation/vcipriani/geneBurdenRD
#BSUB -n 2 
#BSUB -J geneBurdenRD_fisher_jobarray[1-226]
#BSUB -R rusage[mem=10000]
#BSUB -M 10000

module load lang/R/4.1.0-foss-2019b

# Specify 2 mandatory R arguments in quotes:
# args[1]. job index: e.g. 'LSB_JOBINDEX'
# arg[2]. path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"

Rscript scripts/geneBurdenRD_fisher_jobarray.R 'LSB_JOBINDEX' "data/analysisLabelList.tsv"

# To run this script:

# 1. Uncomment: runMode <- "cluster" in geneBurdenRD_fisher_jobarray.R
# 2. ssh cluster
# 3. bsub < /re_gecip/enhanced_interpretation/vcipriani/geneBurdenRD/scripts/geneBurdenRD_fisher_jobarray.sh 
# 4. bjobs 
