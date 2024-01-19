# Specify 1 mandatory R argument:
# args[1]. project name (e.g. same name of project folder) in quotes: "geneBurdenRD"

Rscript scripts/geneBurdenRD_padjust.R "geneBurdenRD"

# To run this script locally, from 'geneBurdenRD' folder:
# 'local' folder has been created before running scripts/geneBurdenRD_prepare_local.sh
# sh scripts/geneBurdenRD_padjust_local.sh > local/geneBurdenRD_padjust.stdout 2> local/geneBurdenRD_padjust.stderr
