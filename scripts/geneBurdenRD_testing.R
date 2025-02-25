
# START ====

args <- commandArgs(trailingOnly = TRUE)

# Load libraries ---- 
library("tidyverse")
library("logistf")

# Get the jobindex ----
jobindex <- as.numeric(args[1]) # <-- user-provided jobindex (e.g. LSB_JOBINDEX)

plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

output <- "output"
dir.create(file.path("./", output), showWarnings = FALSE)

testing <- "testing"
dir.create(file.path("./output/", testing), showWarnings = FALSE)

analysisLabelListFile <- args[2]  # e.g. "data/analysisLabelList.tsv" <--- user-provided file

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
  
  if (length(args) == 5) {
    
    covariatesFile <- args[4]  # e.g. "data/covariates.tsv" <--- user-provided file
    
    # Read in covariates ----
    covariates <- read_tsv(covariatesFile)
    covariateList <- args[5]
    covariateList <- str_split(covariateList, ",")
    covariates <- covariates %>%
      select(sample.id, unlist(covariateList))
    
    matrix <- inner_join(matrix, covariates, by = "sample.id")
    
  }
  if (length(args) == 6) {
    
    covariatesFile <- args[4]  # e.g. "data/covariates.tsv" <--- user-provided file
    
    # Read in covariates ----
    covariates <- read_tsv(covariatesFile)
    covariateList <- args[6]
    covariateList <- str_split(covariateList, ",")
    covariates <- covariates %>%
      select(sample.id, unlist(covariateList))
    
    matrix <- inner_join(matrix, covariates, by = "sample.id")

  }

  # order variable names and create startcol to test genes 
  if ("caco.denovo" %in% colnames(matrix)) {
    
    if (length(args) == 5) {
      print("Missing mandatory args[6]")
    }
    
    # adjusted analysis 
    if (length(args) == 6) {
      
      covariateList <- args[6]
      covariateList <- str_split(covariateList, ",")
      
      matrix <- matrix %>%
        select(caco, caco.denovo, sample.id, unlist(covariateList), everything())
      startcol <- 3 + ncol(covariates)
    
    }
    # unadjusted analysis 
    if (length(args) == 3) {
      matrix <- matrix %>%
        select(caco, caco.denovo, sample.id, everything())
      startcol <- 4
    }
    
  } else {

    # adjusted analysis 
    if (length(args) == 5) {
      
      covariateList <- args[5]
      covariateList <- str_split(covariateList, ",")
      
      matrix <- matrix %>%
        select(caco, sample.id, unlist(covariateList), everything())
      startcol <- 2 + ncol(covariates)
    }
    # unadjusted analysis 
    if (length(args) == 3) {
      matrix <- matrix %>%
        select(caco, sample.id, everything())
      startcol <- 3
    }
  }
  
  # Create empty datalist and info dataframe to populate with loop per column
  datalist = list()
  info <- as.data.frame(matrix(ncol = 0, nrow = 1))
  
  # Loop per each gene and each test ---- 
  if (sum(matrix$caco, na.rm = TRUE) >= 5) { # Test only if total number disease cases is >= 5 
  
    for (k in colnames(matrix)[startcol:ncol(matrix)]){

      colnames_denovo <- names(matrix %>% select(contains("denovo")))
      if (k %in% colnames_denovo == TRUE) {
        
        # adjusted analysis 
        if (length(args) == 6) {
          covariateList <- args[6]
          covariateList <- str_split(covariateList, ",")
        
          # Create denovo contingency table 
          matrix_denovo <- matrix %>% select(k, caco.denovo, unlist(covariateList))
          cont_denovo <- as.matrix(table(matrix_denovo[[k]], matrix_denovo$caco.denovo))
        }
        
        # unadjusted analysis 
        if (length(args) == 3) {
          
         # Create denovo contingency table 
         matrix_denovo <- matrix %>% select(k, caco.denovo)
         cont_denovo <- as.matrix(table(matrix_denovo[[k]], matrix_denovo$caco.denovo))
        }
        
        print(k)
        print("Print caco.denovo contingency table:")
        print(cont_denovo)
        
        if (dim(cont_denovo)[1] == 2) { # Test only if 2x2 contingency table exists
        
        # Criteria to satisfy: at least one disease case with event (d >= 1) AND at least total of 5 disease cases AND at least total of 4 events ----
          if ((cont_denovo[2,2] >= 1 & sum(cont_denovo[,2]) >= 5 & sum(cont_denovo[2,]) >= 4)) {
            
            print("Testing:")
            
            info$analysis.label <- analysis
            info$gene.test <- k
            
            info$totcases <- sum(cont_denovo[,2])
            info$totcontrols <- sum(cont_denovo[,1])
            
            info$a <- cont_denovo[1,1] 
            info$b <- cont_denovo[1,2]
            info$c <- cont_denovo[2,1]
            info$d <- cont_denovo[2,2]
            
            testvar <- args[3]
            
            if (testvar == "fisher") {
              
             # Fisher ----
            
             # pvalue
             info$pvalue <- fisher.test(cont_denovo, alternative = "greater")$p.value
             # OR
             info$or <- (info$a*info$d)/(info$b*info$c)
             # 95% CI
             info$lowcior <- exp(log(info$or) - 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
             info$upcior <-  exp(log(info$or) + 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            
             print("Fisher's test done")
             
            } else if (testvar == "firth") {
              
            # Firth's logistic regression - unadjusted ----
            
            formula <- as.formula(paste0("caco.denovo ~ ", "`", k, "`"))
            
            logistf <- logistf(formula, data = matrix_denovo, control = logistf.control(maxit=25, maxstep=10))
            
            # pvalue
            info$pvalue <- logistf[["prob"]][2]
            # OR
            info$or <- exp(logistf[["coefficients"]][2])
            # 95% CI
            info$lowcior <- exp(logistf[["ci.lower"]][2])
            info$upcior <- exp(logistf[["ci.upper"]][2])
            
            print("Unadjusted Firth's logistic regression done")
            
            } else if (testvar == "adjfirth") {
            
            # Firth's logistic regression - adjusted ----
            
            covariateList <- args[6]
            covariateList <- gsub(",", " + ", covariateList)
            
            formula <- as.formula(paste0("caco.denovo ~ ", "`", k, "`", " + ", unlist(covariateList)))
            print(formula)
            
            logistf <- logistf(formula, data = matrix_denovo, control = logistf.control(maxit=25, maxstep=10))
            
            # pvalue 
            info$pvalue <- logistf[["prob"]][2]
            # OR
            info$or <- exp(logistf[["coefficients"]][2])
            # 95% CI
            info$lowcior <- exp(logistf[["ci.lower"]][2])
            info$upcior <- exp(logistf[["ci.upper"]][2])
            
            print("Firth's adjusted logistic regression done")
            
            }
            else {
              print("Missing mandatory args[3]: fisher/firth/adjfirth")
              
            }
            
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
        
        # adjusted analysis 
        if (length(args) == 5) { 
         
         covariateList <- args[5]
         covariateList <- str_split(covariateList, ",")

         # Create contingency table for tests other than denovo ----
         matrix_not_denovo <- matrix %>% select(k, caco, unlist(covariateList))
         cont <- as.matrix(table(matrix_not_denovo[[k]], matrix_not_denovo$caco))
        }
        if (length(args) == 6) { 
          
          matrix <- read_tsv(analysismatrixFile)
          covariateList <- args[5]
          covariateList <- str_split(covariateList, ",")
          
          covariates <- read_tsv(covariatesFile)
          covariates <- covariates %>%
           select(sample.id, unlist(covariateList))
          matrix <- inner_join(matrix, covariates, by = "sample.id")
          
          # Create contingency table for tests other than denovo ----
          matrix_not_denovo <- matrix %>% select(k, caco, unlist(covariateList))
          cont <- as.matrix(table(matrix_not_denovo[[k]], matrix_not_denovo$caco))
        }
        
        # unadjusted analysis 
        if (length(args) == 3) {
         # Create contingency table for tests other than denovo ----
         matrix_not_denovo <- matrix %>% select(k, caco)
         cont <- as.matrix(table(matrix_not_denovo[[k]], matrix_not_denovo$caco))
        }
        
        print(k)
        print("Print caco contingency table:")
        print(cont)
        
        if (dim(cont)[1] == 2) { # Test only if 2x2 contingency table exists 
          
          # Criteria to satisfy: at least one disease case with event (d >= 1) AND at least total of 5 disease cases AND at least total of 4 events ----
          if ((cont[2,2] >= 1 & sum(cont[,2]) >= 5 & sum(cont[2,]) >= 4)) {
            
            print("Testing:")

            info$analysis.label <- analysis
            info$gene.test <- k
            
            info$totcases <- sum(cont[,2])
            info$totcontrols <- sum(cont[,1])
            
            info$a <- cont[1,1] 
            info$b <- cont[1,2]
            info$c <- cont[2,1]
            info$d <- cont[2,2]
            
            
            testvar <- args[3]
            
            if (testvar == "fisher") {
              
              # Fisher ---- 
            
              # pvalue
              info$pvalue <- fisher.test(cont, alternative = "greater")$p.value
              # OR
              info$or <- (info$a*info$d)/(info$b*info$c)
              # 95% CI
              info$lowcior <- exp(log(info$or) - 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
              info$upcior <-  exp(log(info$or) + 1.96 * sqrt((1/info$a + 1/info$b + 1/info$c + 1/info$d)))
            
              print("Fisher's test done")
            
            } else if (testvar == "firth") { 
              
            # Firth's logistic regression - unadjusted ----
            
            formula <- as.formula(paste0("caco ~ ", "`", k, "`"))
            
            logistf <- logistf(formula, data = matrix_not_denovo, control = logistf.control(maxit=25, maxstep=10))
            
            # pvalue
            info$pvalue <- logistf[["prob"]][2]
            # OR
            info$or <- exp(logistf[["coefficients"]][2])
            # 95% CI
            info$lowcior <- exp(logistf[["ci.lower"]][2])
            info$upcior <- exp(logistf[["ci.upper"]][2])
            
            print("Unadjusted Firth's logistic regression done")
            
            
            } else if (testvar == "adjfirth") {
            
            # Firth's logistic regression - adjusted ----

            covariateList <- args[5]
            covariateList <- gsub(",", " + ", covariateList)
            
            formula <- as.formula(paste0("caco ~ ", "`", k, "`", " + ", unlist(covariateList)))
            print(formula)
            
            #logistf <- logistf(formula, data = matrix_not_denovo)
            logistf <- logistf(formula, data = matrix_not_denovo, control = logistf.control(maxit=30, maxstep=10))
            info$pvalue <- logistf[["prob"]][2]
            info$or <- exp(logistf[["coefficients"]][2])
            info$lowcior <- exp(logistf[["ci.lower"]][2])
            info$upcior <- exp(logistf[["ci.upper"]][2])

            print("Adjusted Firth's logistic regression done")
            }
            else {
              print("Missing mandatory args[3]: fisher/firth/adjfirth")
              
            }
            
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
    outputresults <- paste0("output/testing/", analysis, ".tsv") 
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
