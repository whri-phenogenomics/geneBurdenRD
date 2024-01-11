
# Use R/4.1.0

# START ====

args <- commandArgs(trailingOnly = TRUE)

# Select runMode, either local-interactive or cluster ----
#runMode <- "local-inter"
runMode <- "cluster"

# Load libraries ---- 
library("tidyverse")
library("data.table")
library("reshape2")

if (runMode == "local-inter") {
  
  # Define the jobindex as # of analysis interactively (for example, jobindex <- 1 to select the first analysis) ----
  jobindex <- 1
  
} else {

  # Get the jobindex from cluster ----
  # jobindex <- as.numeric(Sys.getenv('LSB_JOBINDEX')) # <-- user-provided jobindex (e.g. LSB_JOBINDEX)
  jobindex <- as.numeric(Sys.getenv(args[1])) # <-- user-provided jobindex (e.g. LSB_JOBINDEX)
  
}

plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

output <- "output"
dir.create(file.path("./", output), showWarnings = FALSE)

matrix <- "matrix"
dir.create(file.path("./output/", matrix), showWarnings = FALSE)

cluster <- "cluster"
dir.create(file.path("./", cluster), showWarnings = FALSE)



if (runMode == "local-inter") {
  
  exomiserFile <- "data/exomiserPassWide.tsv"
  
} else {
  
  exomiserFile <- args[2]  # e.g. "data/exomiserPassWide.tsv" <--- either user-provided file, or output file from geneBurdenRD_prepare.R/.sh scripts
  
}

if (runMode == "local-inter") {
  
  analysisLabelListFile <- "data/analysisLabelList.tsv"
  
} else {
  
  analysisLabelListFile <- args[3]  # e.g. "data/analysisLabelList.tsv" <--- user-provided file
  # analysisLabelList.tsv is a user-provided file: no-header, 1 column, containing no-spaced analysis labels (e.g. CVD) per each case-control analysis.
  # This file is accompanied by corresponding user-provided case-control files named diseaseLabel.tsv (e.g. CVD.tsv) containing two columns as:
  # sample.id and caco (0/1/NA) where 0 is control, 1 is case, NA missing disease status.
  # The number of sample ids match the number of sample ids run on Exomiser
  
}

# Load exomiserPassWide dataset ---- 
# This file is created in geneBurdenRD_prepare.R 
exomiser <- read_tsv(exomiserFile)

# Create counter variable by sample.id
exomiser <- exomiser %>%
  group_by(sample.id) %>%
  mutate(counter = row_number()) 
# Ungroup 
exomiser <- exomiser %>%
  ungroup()

# Check number of sample ids
length(exomiser$sample.id %>%
         unique()
)

# Get list of analysis labels ----
analyses <- read_tsv(analysisLabelListFile)

# Show number of analyses
dim(analyses)

# If analysisLabelListFile contains only one disease directly impose jobindex as 1
if(nrow(analyses)==1){
  
  jobindex <- 1
  
}

# Define analysis using the jobindex as # of analysis to test  ----
analysis <- analyses[jobindex, ] %>% .[[1]]
analysis

# Import caco variable ---- 
cacoFile <- paste0("data/", analysis, ".tsv") 
caco <- read_tsv(cacoFile)
dim(caco)
print("Print caco table")
table(caco$caco)

# merge case-control definition
exomiser <- left_join(exomiser, caco, by = "sample.id")  
colnames(exomiser)

# Check number of sample ids
print("Number of sample IDs")
length(exomiser$sample.id %>%
         unique()
)

# Print caco table(s) 
print("Print caco table")
exomiser %>% 
  group_by(caco) %>% filter(counter == 1) %>%
  summarise(n = n())
if ("caco.denovo" %in% colnames(exomiser)) {
  print("Print caco.denovo table")
  exomiser %>% 
    group_by(caco.denovo) %>% filter(counter == 1) %>%
    summarise(n = n())
}

# Create variable for testing ----
# Create 0|1 for zero80 
exomiser$zero80 <- ifelse(exomiser$max.varscore >= 0.8, 1, 0)
# Create 0|1 for LoF
exomiser$LoF <- ifelse(exomiser$max.vartype == "LoF", 1, 0)
# Create 0|1 for ccr (*with max.varscore >= 0.8)
exomiser$ccr <- ifelse(exomiser$max.ccr == 1 & exomiser$max.varscore >= 0.8, 1, 0)
if ("caco.denovo" %in% colnames(exomiser)) {
# Create 0|1 for denovo (*family structure must be >= 3)
exomiser$denovo <- ifelse(exomiser$max.denovo == 1 & exomiser$fam.structure >= 3, 1, 0)
}

# Create list of genes where there is at least one proband for the test disease *with the event of interest 
# zero80 
genes80 <- exomiser %>% 
  filter(caco == 1 & max.varscore >= 0.8) %>% 
  select(gene) %>% 
  unique() %>%
  .[[1]]
# LoF 
genesLoF <- exomiser %>% 
  filter(caco == 1 & max.vartype == "LoF") %>% 
  select(gene) %>% 
  unique() %>%
  .[[1]]
# ccr 
genesccr <- exomiser %>% 
  filter(caco == 1 & max.ccr == 1 & max.varscore >= 0.8) %>% 
  select(gene) %>%
  unique() %>%
  .[[1]]
# denovo
if ("caco.denovo" %in% colnames(exomiser)) {
  genesdenovo <- exomiser %>% 
    filter(caco.denovo == 1 & max.denovo == 1 & fam.structure >= 3) %>% 
    select(gene) %>%
    unique() %>%
    .[[1]]
} else {
  genesdenovo <- character() # create an empty genelist
}

# Intersect the gene lists with user gene list (if provided as argument) ----
if (length(args) == 4) {
  
  geneListFileON <- args[4]
  print(geneListFileON)
  
  geneListFile <- paste0("data/", analysis, "_genelist", ".tsv") 
  genelist <- read_tsv(geneListFile, col_names = F) %>% 
    .[[1]]
  
  if (length(genelist) != 0) {
    
    genes80 <- intersect(genes80, genelist)
    genesLoF <- intersect(genesLoF, genelist)
    
    if("caco.denovo" %in% colnames(exomiser)) {
      genesdenovo <- intersect(genesdenovo, genelist)
    }
    
    genesccr <- intersect(genesccr, genelist)
  }
}

# Create zero80 matrix ####
if (length(genes80) >= 1) {
  zero80 <- reshape2::dcast(exomiser, sample.id ~ gene, value.var = 'zero80', fill = 0, drop = FALSE)
  zero80 <- zero80 %>% select(sample.id, genes80)
  colnames(zero80)[2:ncol(zero80)] <- paste(colnames(zero80)[2:ncol(zero80)], "zero80", sep = "_")
  print("Check NA in zero80")
  sum(is.na(zero80)) 
}
# Create LoF  matrix ####
if (length(genesLoF) >= 1) {
  LoF <- reshape2::dcast(exomiser, sample.id ~ gene, value.var = 'LoF', fill = 0, drop = FALSE)
  LoF <- LoF %>% select(sample.id, genesLoF)
  colnames(LoF)[2:ncol(LoF)] <- paste(colnames(LoF)[2:ncol(LoF)], "LoF", sep = "_")
  print("Check NA in LoF")
  sum(is.na(LoF))
}
# Create ccr matrix  ####
if (length(genesccr) >= 1) {
  ccr <- reshape2::dcast(exomiser, sample.id ~ gene, value.var = 'ccr', fill = 0, drop = FALSE)
  ccr <- ccr %>% select(sample.id, genesccr)
  colnames(ccr)[2:ncol(ccr)] <- paste(colnames(ccr)[2:ncol(ccr)], "ccr", sep = "_")
  print("Check NA in ccr")
  sum(is.na(ccr))
} 
# Create denovo matrix  ####
if (length(genesdenovo) >= 1) {
  denovo <- reshape2::dcast(exomiser, sample.id ~ gene, value.var = 'denovo', fill = 0, drop = FALSE)
  denovo <- denovo %>% select(sample.id, genesdenovo)
  colnames(denovo)[2:ncol(denovo)] <- paste(colnames(denovo)[2:ncol(denovo)], "denovo", sep = "_")
  print("Check NA in denovo")
  sum(is.na(denovo))
}

# Merge all gene/test matrixes and create final matrix ---- 
if (length(genes80) >= 1) {
  
  if (length(genesLoF) >= 1) {
    
    matrix <- merge(x = zero80, y = LoF, by = "sample.id", all.x = TRUE, all.y = TRUE)
    
    if (length(genesdenovo) >= 1) {
      
      matrix <- merge(x = matrix, y = denovo, by = "sample.id", all.x = TRUE, all.y = TRUE)
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = matrix, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
      }
      
    } else {
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = matrix, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
      }
    }
    
  } else {
    
    if (length(genesdenovo) >= 1) {
      
      matrix <- merge(x = zero80, y = denovo, by = "sample.id", all.x = TRUE, all.y = TRUE)
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = matrix, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
      }
      
    } else {
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = zero80, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
        
      } else {
        
        matrix <- zero80 
      }
    }
  }
  
  # Get unique proband.info: sample.id and caco, caco.denovo (if exists)
  # This can be written differently in case 
  proband.info <- exomiser[!duplicated(exomiser[ , c('sample.id')]), ]
  matrix$caco <- proband.info$caco[match(matrix$sample.id, proband.info$sample.id)]
  
  if("caco.denovo" %in% colnames(exomiser)) {
    
    matrix$caco.denovo <- proband.info$caco.denovo[match(matrix$sample.id, proband.info$sample.id)]
    matrix <- matrix %>% select(caco, caco.denovo, everything())
    
  } else {
    
    matrix <- matrix %>% select(caco, everything())
  }
    
  # Output disease matrix file as .tsv.gz ----
  # Needs to be compressed to save space 
  outputmatrix <- paste0("output/matrix/", analysis, ".tsv.gz") 
  write_tsv(matrix, gzfile(outputmatrix))
  print("Output disease matrix file as .tsv.gz: completed")
  
} else {
  
  if (length(genesLoF) >= 1) {
    
    if (length(genesdenovo) >= 1) {
      matrix <- merge(x = LoF, y = denovo, by = "sample.id", all.x = TRUE, all.y = TRUE)
      
      if (length(genesccr) >= 1) {
        matrix <- merge(x = matrix, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
      }
    } else {
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = LoF, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
        
      } else {
        
        matrix <- LoF
      }
    }
    
    # Get unique proband.info: sample.id and caco, caco.denovo 
    # This can be written differently in case 
    proband.info <- exomiser[!duplicated(exomiser[ , c('sample.id')]), ]
    
    matrix$caco <- proband.info$caco[match(matrix$sample.id, proband.info$sample.id)]

    if("caco.denovo" %in% colnames(exomiser)) {
      
      matrix$caco.denovo <- proband.info$caco.denovo[match(matrix$sample.id, proband.info$sample.id)]
      matrix <- matrix %>% select(caco, caco.denovo, everything())
      
    } else {
      
      matrix <- matrix %>% select(caco, everything())
      
    }
    # Output disease matrix file as .tsv.gz ----
    # Needs to be compressed to save space 
    outputmatrix <- paste0("output/matrix/", analysis, ".tsv.gz") 
    write_tsv(matrix, gzfile(outputmatrix))
    print("Output disease matrix file as .tsv.gz: completed")
    
  } else {
    
    if (length(genesdenovo) >= 1) {
      
      if (length(genesccr) >= 1) {
        
        matrix <- merge(x = denovo, y = ccr, by = "sample.id", all.x = TRUE, all.y = TRUE)
        
      } else {
        
        matrix <- denovo 
      }
      
      # Get unique proband.info: sample.id and caco, caco.denovo 
      # This can be written differently in case 
      proband.info <- exomiser[!duplicated(exomiser[ , c('sample.id')]), ]
      
      matrix$caco <- proband.info$caco[match(matrix$sample.id, proband.info$sample.id)]
      
      if("caco.denovo" %in% colnames(exomiser)) {
        
        matrix$caco.denovo <- proband.info$caco.denovo[match(matrix$sample.id, proband.info$sample.id)]
        matrix <- matrix %>% select(caco, caco.denovo, everything())
      } else {
        
        matrix <- matrix %>% select(caco, everything())
        
      }
      # Output disease matrix file as .tsv.gz ----
      # Needs to be compressed to save space 
      outputmatrix <- paste0("output/matrix/", analysis, ".tsv.gz") 
      write_tsv(matrix, gzfile(outputmatrix))
      print("Output disease matrix file as .tsv.gz: completed")
      
    } else {
      
      if (length(genesccr) == 0) {
        
        print("Print genes80")
        print(genes80)
        print("Print genesLoF")
        print(genesLoF)
        print("Print genesdenovo")
        print(genesdenovo)
        print("Print genesccr")
        print(genesccr)
        
        print("No output disease matrix file")
        
      } else {
        
        matrix <- ccr 
        
        # Get unique proband.info: sample.id and caco, caco.denovo 
        # This can be written differently in case 
        proband.info <- exomiser[!duplicated(exomiser[ , c('sample.id')]), ]
        
        matrix$caco <- proband.info$caco[match(matrix$sample.id, proband.info$sample.id)]
        
        if("caco.denovo" %in% colnames(exomiser)) {
          
          matrix$caco.denovo <- proband.info$caco.denovo[match(matrix$sample.id, proband.info$sample.id)]
          matrix <- matrix %>% select(caco, caco.denovo, everything())
          
        } else {
          
          matrix <- matrix %>% select(caco, everything())
          
        } 
        # Output disease matrix file as .tsv.gz ----
        # Needs to be compressed to save space 
        outputmatrix <- paste0("output/matrix/", analysis, ".tsv.gz") 
        write_tsv(matrix, gzfile(outputmatrix))
        print("Output disease matrix file as .tsv.gz: completed")
        
      }
    }
  }
}

# Clean ====
if (runMode == "local-inter") {
  cat("\014")
}

rm(list = ls())

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#           ==== End ==== 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
