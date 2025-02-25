# Specify 2 mandatory R arguments:
# args[1]. loop index that will go from 1 to total number of case-control analyses: e.g. 'I'
# arg[2]. path for file containing list of case-control analysis labels: e.g. "data/analysisLabelList.tsv"
# (path relative to cwd; otherwise, absolute path; in quotes)
# args[3]. specify the type of statistical test among "fisher", "firth", "adjfirth": e.g. "adjfirth"
# "fisher" indicates one-sided Fisher's test; "firth" indicates Firth's logistic regression;
# "adjfirth" indicates covariate-adjusted Firth's logistic regression
# IMPORTANT - if "adjfirth" is selected, 2 more mandatory arguments are needed:
# args[4]. path for file containing sample ids and any covariates: e.g. "data/covariates.tsv"
# args[5]. list of comma-separated (no spaces) covariate variable names consistent with header of covariates file:
# e.g. "age,Participant.Phenotypic.Sex,inferred_anc,family.structure"
# Please note if categorical variables are coded as numbers, specify them as factor(variablename):
# e.g. "age,factor(Participant.Phenotypic.Sex),inferred_anc,factor(family.structure)"
# One last mandatory argument, if de-novo tests are performed (when caco.denovo exists in the user-provided case-control files named analysis.label.tsv)
# and "adjfirth" is selected, provide:
# args[6]. list of comma-separated (no spaces) covariate variable names consistent with header of covariates file;
# this is either the same as above or a different, specific list of covariates for de-novo tests only
# e.g. "age,Participant.Phenotypic.Sex,inferred_anc"
# Please note if categorical variables are coded as numbers, specify them as factor(variablename):
# e.g. "age,factor(Participant.Phenotypic.Sex),inferred_anc"

# Loop over all case-control analyses in "data/analysisLabelList.tsv"

TOTAL=$1;
if [ -z "$1" ]
then
        TOTAL=1;
fi
for I in `seq 1 $TOTAL`; do
   	echo ${I}
   	Rscript scripts/geneBurdenRD_testing.R ${I} "data/analysisLabelList.tsv" "adjfirth" "data/sample_covariates.tsv" "age,Participant.Phenotypic.Sex,inferred_anc,family.structure" "age,Participant.Phenotypic.Sex,inferred_anc"
	#Rscript scripts/geneBurdenRD_testing.R ${I} "data/analysisLabelList.tsv" "fisher"	
        #Rscript scripts/geneBurdenRD_testing.R ${I} "data/analysisLabelList.tsv" "firth"
	
done

# To run this script locally, from 'geneBurdenRD' folder:
# 'local' folder has been created before running scripts/geneBurdenRD_prepare_local.sh
# sh scripts/geneBurdenRD_testing_local.sh > local/geneBurdenRD_testing.stdout 2> local/geneBurdenRD_testing.stderr
