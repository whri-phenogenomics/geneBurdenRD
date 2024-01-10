#BSUB -q long
#BSUB -P re_gecip_enhanced_interpretation
#BSUB -o cluster/geneBurdenRD_visualisation_jobarray_%J_%I.stdout
#BSUB -e cluster/geneBurdenRD_visualisation_jobarray_%J_%I.stderr
#BSUB -cwd /geneBurdenRD
#BSUB -n 5
#BSUB -J geneBurdenRD_visualisation_jobarray[1-226]
#BSUB -R rusage[mem=300000]
#BSUB -M 300000

module load lang/R/4.1.0-foss-2019b

# Specify 3 mandatory R arguments in quotes:
# args[1]. job index: e.g. 'LSB_JOBINDEX'
# args[2]. path for Exomiser masterfile file: e.g. "data/exomiser_master_file_passvars.tsv"
# args[3]. path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"

# One optional argument
# args[4]. choice of FDR threshold: "0.01" default is "0.05"
# args[5]. choice "plot" or "table_only" default is "plot"

Rscript /geneBurdenRD/scripts/geneBurdenRD_visualisation_jobarray.R 'LSB_JOBINDEX' "data/exomiser_master_file_passvars.tsv" "data/analysisLabelList.tsv" "0.05" "plot"
