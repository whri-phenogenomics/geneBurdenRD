
# Use R/4.1.0

# START ====

args <- commandArgs(trailingOnly = TRUE)

# Load libraries ---- 
library("tidyverse")

# Get the jobindex ----
jobindex <- as.numeric(args[1]) # <-- user-provided jobindex (e.g. LSB_JOBINDEX)

plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

output <- "output"
dir.create(file.path("./", output), showWarnings = FALSE)

fisher <- "fisher"
dir.create(file.path("./output/", fisher), showWarnings = FALSE)

cluster <- "cluster"
dir.create(file.path("./", cluster), showWarnings = FALSE)

analysisLabelListFile <- args[2]  # e.g. "data/analysisLabelList.tsv" <--- user-provided file
  # analysisLabelList.tsv is a user-provided file: no-header, 1 column, containing no-spaced analysis labels (e.g. CVD) per each case-control analysis.
  # This file is accompanied by corresponding user-provided case-control files named analysis.label.tsv (e.g. CVD.tsv) containing two columns as:
  # sample.id and caco (0/1/NA) where 0 is control, 1 is case, NA missing disease status.
  # The number of sample ids match the number of sample ids run on Exomiser

# Get list of analysis labels ----
analyses <- read_tsv(analysisLabelListFile)

# Show number of analyses
dim(analyses)

# Define analysis using the jobindex as # of analysis to test  ----
analysis <- analyses[jobindex, ] %>% .[[1]]
analysis

analysisLabel <- analyses[jobindex, ] %>% .[[2]]
analysisLabel

# Read in analysis matrix ----
analysismatrixFile <- paste("output/matrix/", analysis, ".tsv.gz", sep = "") 

if (file.exists(analysismatrixFile) == TRUE ) {
  
  matrix <- read_tsv(analysismatrixFile)
  str(matrix)

  # Create empty datalist and info dataframe to populate with loop per column
  datalist = list()
  info <- as.data.frame(matrix(ncol = 0, nrow = 1))
  
  # Loop per each gene and each test ---- 
  if (sum(matrix$caco, na.rm = TRUE) >= 5) { # Test only if total number disease cases is >= 5 
  
    for (k in colnames(matrix)[4:ncol(matrix)]){

      colnames_denovo <- names(matrix %>% select(contains("denovo")))
      if (k %in% colnames_denovo == TRUE) {
        
        # Create denovo contingency table 
        matrix_denovo <- matrix %>% select(k, caco.denovo)
        cont_denovo <- as.matrix(table(matrix_denovo[[k]], matrix_denovo$caco.denovo))
        print(cont_denovo)
        
        if (dim(cont_denovo)[1] == 2) { # Test only if 2x2 contingency table exists
        
        # Criteria to satisfy: at least one disease case with event (d >= 1) AND at least total of 5 disease cases AND at least total of 4 events ----
          if ((cont_denovo[2,2] >= 1 & sum(cont_denovo[,2]) >= 5 & sum(cont_denovo[2,]) >= 4)) {
            
            print("Testing:")
            print(k)
            
            info$analysis.label <- analysis
            info$gene.test <- k
            
            info$pvalue <- fisher.test(cont_denovo, alternative = "greater")$p.value

            info$totcases <- sum(cont_denovo[,2])
            info$totcontrols <- sum(cont_denovo[,1])
            
            info$a <- cont_denovo[1,1] 
            info$b <- cont_denovo[1,2]
            info$c <- cont_denovo[2,1]
            info$d <- cont_denovo[2,2]
            
            # OR
            info$or <- (info$a*info$d)/(info$b*info$c)
            # 95% CI
            info$lowcior <- exp(log(info$or) - 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            info$upcior <-  exp(log(info$or) + 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            
            print("Fisher's test done")
            
            # Fill in datalist 
            datalist[[k]] <- info
    
          } else {
            
            print(k)
            print("2x2 exists, but criteria not satisfied")
            
          }
        } else {
          
          print(k)
          print("Not a 2x2 contingency table")
          
        }
        
      } else {
        
        # Create contingency table for tests other than denovo ----
        matrix_not_denovo <- matrix %>% select(k, caco)
        cont <- as.matrix(table(matrix_not_denovo[[k]], matrix_not_denovo$caco))
        print(cont)
        
        if (dim(cont)[1] == 2) { # Test only if 2x2 contingency table exists 
          
          # Criteria to satisfy: at least one disease case with event (d >= 1) AND at least total of 5 disease cases AND at least total of 4 events ----
          if ((cont[2,2] >= 1 & sum(cont[,2]) >= 5 & sum(cont[2,]) >= 4)) {
            
            print(k)
            
            info$analysis.label <- analysis
            info$gene.test <- k
            
            info$pvalue <- fisher.test(cont, alternative = "greater")$p.value

            info$totcases <- sum(cont[,2])
            info$totcontrols <- sum(cont[,1])
            
            info$a <- cont[1,1] 
            info$b <- cont[1,2]
            info$c <- cont[2,1]
            info$d <- cont[2,2]
            
            # OR
            info$or <- (info$a*info$d)/(info$b*info$c)
            # 95% CI
            info$lowcior <- exp(log(info$or) - 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            info$upcior <-  exp(log(info$or) + 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            
            print("Fisher's test done")
            
            # Fill in datalist 
            datalist[[k]] <- info
    
          } else {
            
            print(k)
            print("2x2 exists, but criteria not satisfied")
            
          }
        } else {
          
          print(k)
          print("Not a 2x2 contingency table")
          
        }
    
      }
    }
    
    # Collect all the results for the disease tested ----
    results = do.call(rbind, datalist)
    results$analysis.label <- analysis
    results$analysis <- analysisLabel
    results$split <- results$gene.test
    results$tottests <- nrow(results)
    
    results <- separate(results, split, c("gene", "test"), sep = '_') 
    
    results <- results %>% select(analysis.label, analysis, gene, test, pvalue, or, lowcior, upcior, d, totcases, totcontrols, tottests, a, b, c)
    options(scipen = 999)
    
    # Output all the results for the disease tested ----
    outputresults <- paste0("output/fisher/", analysis, ".tsv") 
    write_tsv(results, outputresults)
    print("Output test result file as .tsv: completed") 
  
  } else {
    
    print("Total number of cases is less than 5: disease NOT tested")
  }
  
} else {
  
  print(analysismatrixFile)
  print("This matrix file does NOT exist")
}

rm(list = ls())

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#           ==== End ==== 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
