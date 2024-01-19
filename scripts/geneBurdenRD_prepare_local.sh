# Specify 2 mandatory R argument:
# arg[1]. path for Exomiser/Genomiser candidate variant file (all cases and controls): e.g. "data/exomiser_master_file_passvars.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)
# arg[2]. path for final output file which will be the input data file for matrix script: e.g. "data/exomiserPassWide.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)

Rscript scripts/geneBurdenRD_prepare.R "data/exomiser_master_file_passvars.tsv" "data/exomiserPassWide.tsv"

# To run this script locally, from 'geneBurdenRD' folder:
# Create 'local' folder to redirect the standard output and standard error
# mkdir local
# sh scripts/geneBurdenRD_prepare_local.sh > local/geneBurdenRD_prepare.stdout 2> local/geneBurdenRD_prepare.stderr
