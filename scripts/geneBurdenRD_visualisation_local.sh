# Specify 3 mandatory R arguments:
# args[1]. loop index that will go from 1 to total number of case-control analyses: e.g. 'I'
# args[2]. path for Exomiser masterfile file: e.g. "data/exomiser_master_file_passvars.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)
# args[3]. path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)

# Two optional arguments
# args[4]. choice of FDR threshold: "0.01" - default is "0.05"
# args[5]. choice of "plot" or "table_only" - default is "plot"

# Loop over all case-control analyses in "data/analysisLabelList.tsv"

TOTAL=$1;
if [ -z "$1" ]
then
        TOTAL=1;
fi
for I in `seq 1 $TOTAL`; do
   echo ${I}
   Rscript scripts/geneBurdenRD_visualisation.R ${I} "data/exomiser_master_file_passvars.tsv" "data/analysisLabelList.tsv" "0.05" "plot"

done

# To run this script locally, from 'geneBurdenRD' folder:
# 'local' folder has been created before running scripts/geneBurdenRD_prepare_local.sh
# sh scripts/geneBurdenRD_visualisation_local.sh > local/geneBurdenRD_visualisation.stdout 2> local/geneBurdenRD_visualisation.stderr
