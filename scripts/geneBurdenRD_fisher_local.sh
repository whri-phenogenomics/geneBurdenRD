# Specify 2 mandatory R arguments:
# args[1]. loop index that will go from 1 to total number of case-control analyses: e.g. 'I'
# arg[2]. path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)

# Loop over all case-control analyses in "data/analysisLabelList.tsv"
for I in (1..1}; do

   echo ${I}
   Rscript scripts/geneBurdenRD_fisher.R ${I} "data/analysisLabelList.tsv"

done

# To run this script locally, from 'geneBurdenRD' folder:
# 'local' folder has been created before running scripts/geneBurdenRD_prepare_local.sh
# sh scripts/geneBurdenRD_fisher_local.sh > local/geneBurdenRD_fisher.stdout 2> local/geneBurdenRD_fisher.stderr
