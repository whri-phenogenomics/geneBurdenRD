#BSUB -q short
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -o cluster/geneBurdenRD_padjust_%J.stdout
#BSUB -e cluster/geneBurdenRD_padjust_%J.stderr
#BSUB -cwd /re_gecip/enhanced_interpretation/vcipriani/geneBurdenRD
#BSUB -n 1
#BSUB -R rusage[mem=20000]
#BSUB -M 20000

module load lang/R/4.1.0-foss-2019b

# Specify 1 mandatory R argument in quotes:
# args[1]. project name (e.g. same name of project folder): "geneBurdenRD"

Rscript scripts/geneBurdenRD_padjust.R "geneBurdenRD"

# To run this script:

# 1. ssh cluster
# 2. bsub < /re_gecip/enhanced_interpretation/vcipriani/geneBurdenRD/scripts/geneBurdenRD_padjust.sh 
# 3. bjobs 
