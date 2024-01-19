# Specify 3 mandatory R arguments:
# args[1]. loop index that will go from 1 to total number of case-control analyses: e.g. 'I'
# args[2]. data path for Exomiser/Genomiser file: e.g. "data/exomiserPassWide.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)
# args[3]. data path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)

# One optional argument
# args[4]. whether gene list files are provided: "geneListFileON"

# Loop over all case-control analyses in "data/analysisLabelList.tsv"
for I in (1..1}; do

   echo ${I}
   Rscript scripts/geneBurdenRD_matrix.R ${I} "data/exomiserPassWide.tsv" "data/analysisLabelList.tsv"

done

# To run this script locally, from 'geneBurdenRD' folder:
# 'local' folder has been created before running scripts/geneBurdenRD_prepare_local.sh 
# sh scripts/geneBurdenRD_matrix_local.sh > local/geneBurdenRD_matrix.stdout 2> local/geneBurdenRD_matrix.stderr
