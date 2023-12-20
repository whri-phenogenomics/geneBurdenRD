
# Use R/4.1.0

# START ====

# Load libraries ---- 
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

cluster <- "cluster"
dir.create(file.path("./", cluster), showWarnings = FALSE)

results <- "results"
dir.create(file.path("./", results), showWarnings = FALSE)

# Create read function ----
read_fun <- function(flnm) {
read_tsv(flnm, col_types = c("ccccddddddddddd"), col_names = TRUE) %>% 
  mutate(filename = flnm)
}

# Get list of p-value disease files ----
pvalFileList <- list.files(path = "output/fisher", pattern = "*.tsv", full.names = T)
# Display how many diseases were tested
print("Total number of diseases tested:")
length(pvalFileList)

# Read in all p-value disease files and merge them in a tibble ----
pvalues <- pvalFileList %>%
  map_df(~ read_fun(.))

# Show number of tests 
print("Total number of tests:")
nrow(pvalues)

# FDR adjustment of all p-values (all diseases, all tests) ----
pvalues$p.adjust.fdr <- p.adjust(pvalues$pvalue, "fdr")

# Significant at FDR 0.1 and 0.05
pvalues <- pvalues %>% 
  mutate(signif.fdr.0.10 = ifelse(p.adjust.fdr <= 0.10, 1, 0),
         signif.fdr.0.05 = ifelse(p.adjust.fdr <= 0.05, 1, 0)) 

# Bonferroni correction ----
pvalues <- pvalues %>% 
  mutate(bonferroni.cutoff = 0.05/nrow(pvalues), 
         signif.bonferroni = ifelse(pvalue <= bonferroni.cutoff, 1, 0))

# Total genes tested per each analysis.label
totgenestestedbydisease <-  pvalues %>%
  group_by(analysis.label) %>%
  distinct(gene) %>%
  count()

pvalues <- full_join(pvalues, totgenestestedbydisease, by = "analysis.label")
pvalues <- pvalues %>%
  rename(totgenestested = n)

# Reorder variables  
pvalues <- pvalues %>% select(analysis.label, analysis, gene, test, pvalue, p.adjust.fdr, signif.fdr.0.10, signif.fdr.0.05, or,	d,	totcases,	totcontrols, totgenestested,	tottests, signif.bonferroni, everything())

# Output all adjusted p-value results ---- 
outputadjust <- paste("results/", "geneBurdenRD_FDR_", args[1], ".tsv", sep = "") 
write_tsv(pvalues, outputadjust)
print("Output all adjusted p-value results: completed")

rm(list = ls())

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#           ==== End ==== 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


