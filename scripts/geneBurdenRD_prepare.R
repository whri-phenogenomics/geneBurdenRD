
# Use R/4.1.0

# START ====

args <- commandArgs(trailingOnly = TRUE)

# Load libraries ----
library("tidyverse")
library("data.table")

exomiserCandidatesFile <- args[1] # e.g. "data/exomiser_master_file_passvars.tsv"
  
plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

output <- "output"
dir.create(file.path("./", output), showWarnings = FALSE)

data <- "data"
dir.create(file.path("./", data), showWarnings = FALSE)

# Read in Exomiser candidates file ----
exomiserCandidates <- read_tsv(exomiserCandidatesFile) 

attr(exomiserCandidates, "problems")

str(exomiserCandidates)

# Make variable names and rename variables ====
colnames(exomiserCandidates)
(colnames(exomiserCandidates) <- make.names(tolower(colnames((exomiserCandidates)))))

glimpse(exomiserCandidates)

# Check number of sample ids
# sample.id
length(exomiserCandidates$sample.id %>%
         unique()
)

exomiserPass <- exomiserCandidates
dim(exomiserPass)
rm(exomiserCandidates)

# create counter variable by sample.id
exomiserPass <- exomiserPass %>%
  group_by(sample.id) %>%
  mutate(counter = row_number())
summary(exomiserPass$counter)

exomiserPass <- exomiserPass %>%
  ungroup()

# Check number of sample ids
length(exomiserPass$sample.id %>%
         unique()
)

summary(exomiserPass$variant.score)

# Create unique.casegene.id from sample.id and gene ====
exomiserPass$unique.casegene.id <- paste(exomiserPass$sample.id, exomiserPass$gene, sep = "-")

# Create varcasegene.counter: counts of variants per each sample and each gene ----
# This mutate command may take a while interactively 
exomiserPass <- exomiserPass %>%
  group_by(sample.id, gene) %>%
  mutate(varcasegene.counter = row_number())
exomiserPass <- exomiserPass %>%
  ungroup()  
table(exomiserPass$varcasegene.counter)

# Create vartype from functional.class ==== 
# as presence of:
# at least 1 predicted "LoF" 
# if not, at least 1 "Missense",
# if not, at least 1 "Synonymous",
# if not, "Other" 

exomiserPass$functional.class <- tolower(exomiserPass$functional.class)
exomiserPass$functional.class <- sub(".$", "", exomiserPass$functional.class)

# List functional classes
func.class <- data.frame(table(exomiserPass$functional.class))
func.class <- func.class[order(-func.class$Freq),]

# vartype: at least 1 LoF ----
word <- c("start_lost", "stop_gained", "stop_lost", "frameshift_elongation", "frameshift_truncation", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "exon_loss_variant", "initiator_codon_variant")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiserPass$vartype <- ifelse(grepl(pattern, exomiserPass$functional.class), 'LoF', NA)
func.class$vartype <- ifelse(grepl(pattern, func.class$Var1), 'LoF', NA)

sum(is.na(func.class$vartype))

# vartype: at least 1 missense (if no LoF) ----
word <- c("missense")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiserPass$vartype <- ifelse(is.na(exomiserPass$vartype) == TRUE & grepl(pattern,exomiserPass$functional.class), 'missense', exomiserPass$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern,func.class$Var1), 'missense', func.class$vartype)

sum(is.na(func.class$vartype))

# vartype: at least 1 synonymous (if no LoF, and no missense) ----
word <- c("synonymous")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiserPass$vartype <- ifelse(is.na(exomiserPass$vartype) == TRUE & grepl(pattern,exomiserPass$functional.class), 'synonymous', exomiserPass$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern,func.class$Var1), 'synonymous', func.class$vartype)

sum(is.na(func.class$vartype))

# vartype: other (if no LoF, no missense, no synonymous) ----
word <- c("splice_region", "transcript_ablation", "intron", "upstream" , "sequence", "intergenic", "prime", "downstream", "inframe_insertion", "inframe_deletion", "structural_variant", "stop_retained_variant")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern

exomiserPass$vartype <- ifelse(is.na(exomiserPass$vartype) == TRUE & grepl(pattern, exomiserPass$functional.class), 'other', exomiserPass$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern, func.class$Var1), 'other', func.class$vartype)

sum(is.na(func.class$vartype))

# Check what it is left as NA 
func.class[is.na(func.class$vartype), ]

sum(is.na(exomiserPass$vartype))
func.class %>%
  group_by(vartype) %>%
  summarise(n = n())

rm(func.class)

# Reshape from long to wide ==== 
if ("de.novo" %in% colnames(exomiserPass)) {
  
  wide <- dcast(setDT(exomiserPass),
              unique.casegene.id ~ varcasegene.counter,
              value.var = c("variant.score", "vartype", "de.novo", "ccr.flag", 
                            "variant", "hgvs", "exomisermoi", "genotypes", "functional.class", 
                            "gnomad.freq", "max.freq", "rank"))
  
  test <- exomiserPass[!duplicated(exomiserPass[ , c('unique.casegene.id')]), ]
  dim(test)
  
  # Add other relevant columns (at gene level, or patient level) 
  exomiserPassWide <- merge(x = wide, y = exomiserPass[!duplicated(exomiserPass$unique.casegene.id), c("unique.casegene.id", "sample.id", "gene", "fam.structure")], by = "unique.casegene.id", all.x = TRUE)
  
} else {
  wide <- dcast(setDT(exomiserPass),
                unique.casegene.id ~ varcasegene.counter,
                value.var = c("variant.score", "vartype", "ccr.flag", 
                              "variant", "hgvs", "exomisermoi", "genotypes", "functional.class", 
                              "gnomad.freq", "max.freq", "rank"))
  
  test <- exomiserPass[!duplicated(exomiserPass[ , c('unique.casegene.id')]), ]
  dim(test)
  
  # Add other relevant columns (at gene level, or patient level) 
  exomiserPassWide <- merge(x = wide, y = exomiserPass[!duplicated(exomiserPass$unique.casegene.id), c("unique.casegene.id", "sample.id", "gene")], by = "unique.casegene.id", all.x = TRUE)
  
}

# # Check totals 
table(exomiserPass$varcasegene.counter)
dim(exomiserPassWide) 

# Create new variables in the wide form (per row)

# Create max.denovo ----
if ("de.novo_1" %in% colnames(exomiserPassWide)) {
  
  exomiserPassWide$max.denovo <- ifelse(exomiserPassWide$de.novo_1 == "Y" | exomiserPassWide$de.novo_2 == "Y", "1", "0")
  table(exomiserPassWide$max.denovo)
  exomiserPassWide$max.denovo[is.na(exomiserPassWide$max.denovo)] <- 0
  print("Print max.denovo table")
  table(exomiserPassWide$max.denovo)
  
  exomiserPassWide$de.novo_1 <- NULL
  exomiserPassWide$de.novo_2 <- NULL
}

# Create max.vartype ----
exomiserPassWide$max.vartype <- ifelse(exomiserPassWide$vartype_1 == "LoF" | exomiserPassWide$vartype_2 == "LoF", "LoF", NA)

exomiserPassWide$max.vartype <- ifelse(is.na(exomiserPassWide$max.vartype) == TRUE & (grepl('missense', exomiserPassWide$vartype_1) | grepl('missense', wide$vartype_2)), 'missense', exomiserPassWide$max.vartype)

exomiserPassWide$max.vartype <- ifelse(is.na(exomiserPassWide$max.vartype) == TRUE & (grepl('synonymous', exomiserPassWide$vartype_1) | grepl('synonymous', exomiserPassWide$vartype_2)), 'synonymous', exomiserPassWide$max.vartype)

exomiserPassWide$max.vartype <- ifelse(is.na(exomiserPassWide$max.vartype) == TRUE & (grepl('other', exomiserPassWide$vartype_1) | grepl('other', exomiserPassWide$vartype_2)), 'other', exomiserPassWide$max.vartype)

sum(is.na(exomiserPassWide$max.vartype))
print("Print vartype table")
table(exomiserPassWide$max.vartype)

exomiserPassWide$vartype_1 <- NULL
exomiserPassWide$vartype_2 <- NULL

# Create max.ccr ----
exomiserPassWide$ccr.test_1 <- ifelse((exomiserPassWide$ccr.flag_1 == 1 & exomiserPassWide$variant.score_1 >= 0.80), 1, 0)
exomiserPassWide$ccr.test_2 <- ifelse((exomiserPassWide$ccr.flag_2 == 1 & exomiserPassWide$variant.score_2 >= 0.80), 1, 0)

exomiserPassWide$max.ccr <- pmax(exomiserPassWide$ccr.test_1, exomiserPassWide$ccr.test_2, na.rm = TRUE)
print("Print ccr table")
table(exomiserPassWide$max.ccr)

exomiserPassWide$ccr.flag_1 <- NULL
exomiserPassWide$ccr.flag_2 <- NULL
exomiserPassWide$ccr.test_1 <- NULL
exomiserPassWide$ccr.test_2 <- NULL

# Create max.varscore ---- 
exomiserPassWide$max.varscore <- pmax(exomiserPassWide$variant.score_1, exomiserPassWide$variant.score_2, na.rm = TRUE)

exomiserPassWide$variant.score_1 <- NULL
exomiserPassWide$variant.score_2 <- NULL

# Check
sum(is.na(exomiserPassWide$max.varscore))
if ("max.denovo" %in% colnames(exomiserPassWide)) {
   sum(is.na(exomiserPassWide$max.denovo))
}
sum(is.na(exomiserPassWide$max.vartype))
sum(is.na(exomiserPassWide$max.ccr))

# Select variables to output 
if ("max.denovo" %in% colnames(exomiserPassWide)) {
  
  exomiserPassWideToOutput <- exomiserPassWide %>%
  select(unique.casegene.id, sample.id, gene, fam.structure, max.vartype, max.varscore, max.denovo, max.ccr)
} else {
  exomiserPassWideToOutput <- exomiserPassWide %>%
    select(unique.casegene.id, sample.id, gene,  max.vartype, max.varscore, max.ccr)
}

# Output wide_to_output file ---- 
# This file is used to run the gene-burden testing 
colnames(exomiserPassWideToOutput)
write_tsv(exomiserPassWideToOutput, args[2])
print("Output PassWide.tsv file: completed")

rm(list = ls())

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#           ==== End ==== 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

