#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 
#        Visualisation          #
#      Gene Burdern Test        # 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Use R/4.1.0

args <- commandArgs(trailingOnly = TRUE)

# Load libraries ----
if (! requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
  loadNamespace("tidyverse")
}

library(tidyverse)

#Load results FDR file

path_to_analysis<-getwd()

setwd(path_to_analysis)      # e.g.setwd("/Users/Desktop/geneBurdenRD")

plots <- "plots"
dir.create(file.path("./", plots), showWarnings = FALSE)

file_to_read <-
  list.files(path = "./results",
             pattern = "geneBurdenRD_FDR_.*",
             full.names = T)

w<-str_sub(file_to_read, start=11 ,end=-5)

print(paste0("Reading ",w))

pval_DF<-read_tsv(file_to_read)

#Chose p adjust FDR threshold to use either provided by user or default 0.05

if(length(args)<3){
    
    stop(print(paste("only ",length(args)," arguments provided, please provide [3] mandatory arguments")))

} 

if (length(args) >= 4 && !(args[4] %in% c("plot", "table_only"))) {
  
  print("Argument[4] provided")
  p.adjust.fdr_threshold <-args[4]
  p.adjust.fdr_threshold <-as.numeric(p.adjust.fdr_threshold) #provided by user e.g. "0.01"
  print(paste0("FRD threshold manually set to : ", p.adjust.fdr_threshold))
  
  }else{
  
  p.adjust.fdr_threshold<-as.numeric(0.05)
  print(paste0("FRD threshold set to : ", p.adjust.fdr_threshold," by default"))
  
}


#Discover number of significant gene-disease associations for that threshold
sign_pval_DF<-pval_DF %>% filter(p.adjust.fdr<=p.adjust.fdr_threshold)

#Retrieve additional information from your masterfile
exomiserFile <- args[2] # e.g. "data/exomiser_master_file_passvars.tsv"

exomiser <- read_tsv(exomiserFile)

#Make column name lower case
colnames(exomiser) <- make.names(tolower(colnames(exomiser)))

#Added this only because masterfile contain variants NA
exomiser<-exomiser %>% filter(!is.na(exomiser$variant))

#Clean functional class data
exomiser$functional.class <- tolower(exomiser$functional.class)
exomiser$functional.class <- sub("\\|$", "", exomiser$functional.class)

#Create counter variable by sample.id
exomiser <- exomiser %>%
  group_by(sample.id) %>%
  mutate(counter = row_number())

#Sample list
sampleIDs <- exomiser %>%
  filter(counter == 1) %>%
  dplyr::select(sample.id)

exomiser <- exomiser %>%
  ungroup()

exomiser$unique.casegene.id <- paste(exomiser$sample.id, exomiser$gene, sep = "-")

# Create varcasegene.counter: counts of variants per each GeL patient (case) and each gene
# This mutate command may take a while interactively
exomiser <- exomiser %>%
  group_by(sample.id, gene) %>%
  mutate(varcasegene.counter = row_number()) %>%
  ungroup()

table(exomiser$varcasegene.counter)



# List functional classes
func.class <- data.frame(table(exomiser$functional.class))
func.class <- func.class[order(-func.class$Freq),]

# vartype: at least 1 LoF ----
word <- c("start_lost", "stop_gained", "stop_lost", "frameshift_elongation", "frameshift_truncation", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "exon_loss_variant", "initiator_codon_variant")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiser$vartype <- ifelse(grepl(pattern, exomiser$functional.class), 'LoF', NA)
func.class$vartype <- ifelse(grepl(pattern, func.class$Var1), 'LoF', NA)

sum(is.na(func.class$vartype))

# vartype: at least 1 missense (if no LoF) ----
word <- c("missense")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiser$vartype <- ifelse(is.na(exomiser$vartype) == TRUE & grepl(pattern,exomiser$functional.class), 'missense', exomiser$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern,func.class$Var1), 'missense', func.class$vartype)

sum(is.na(func.class$vartype))

# vartype: at least 1 synonymous (if no LoF, and no missense) ----
word <- c("synonymous")
pattern <- paste0('.*', word, '.*', collapse = '|')
pattern
exomiser$vartype <- ifelse(is.na(exomiser$vartype) == TRUE & grepl(pattern,exomiser$functional.class), 'synonymous', exomiser$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern,func.class$Var1), 'synonymous', func.class$vartype)

sum(is.na(func.class$vartype))

# vartype: other (if no LoF, no missense, no synonymous) ----
#word <- c("splice_region", "transcript_ablation", "intron", "upstream" , "sequence", "internal_feature_elongation", "intergenic", "prime", "downstream", "inframe_insertion", "inframe_deletion", "structural_variant", "stop_retained_variant")
#pattern <- paste0('.*', word, '.*', collapse = '|')
#pattern

#exomiser$vartype <- ifelse(is.na(exomiser$vartype) == TRUE & grepl(pattern, exomiser$functional.class), 'other', exomiser$vartype)
#func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE & grepl(pattern, func.class$Var1), 'other', func.class$vartype)

#To be more flexible I think we should classify anything that is NA into "other"
exomiser$vartype <- ifelse(is.na(exomiser$vartype) == TRUE, 'other', exomiser$vartype)
func.class$vartype <- ifelse(is.na(func.class$vartype) == TRUE, 'other', func.class$vartype)

# Check what it is left as NA 
sum(is.na(func.class$vartype))
sum(is.na(exomiser$vartype))

func.class %>%
  group_by(vartype) %>%
  summarise(n = n())

rm(func.class)

#Filter significant FDR signals

sign_pval_DF<-pval_DF %>% filter(p.adjust.fdr<=p.adjust.fdr_threshold)

print(paste0("There are ",nrow(sign_pval_DF)," significant signals with p.adjust.fdr <= ",p.adjust.fdr_threshold))


#Extract list of unique disease level per significant signals
sign_analysis_label <- setNames(as.data.frame(table(sign_pval_DF$analysis.label)), c("analysis.label", "Freq"))

#Load disease list
analysisLabelListFile <- args[3]  # e.g. "data/analysisLabelList.tsv"
analyses <- read_tsv(analysisLabelListFile)

# Get the jobindex ----
jobindex <- as.numeric(args[1]) # <-- user-provided jobindex (e.g. LSB_JOBINDEX)

disease <- analyses[jobindex, ] %>% .[[1]]

#Stop analysis if disease label has no significant signals
if(!disease %in% sign_analysis_label$analysis.label){
  
  print(paste0("no significant signals for disease label ",disease, " using FDR threshold <= ",p.adjust.fdr_threshold))
  
  break

}


#if plotting is requested as argument -

if (length(args) == 5 && (args[5] %in% c("plot", "table_only"))) {
  
  plotting_arg<-args[5]
  plotting_arg<-as.character(plotting_arg)
  plotting<-ifelse(plotting_arg=="table_only","N","Y")
  
  if((plotting_arg=="table_only")==TRUE){
    
    print("Data will be generated only as table")
    
  }

    
  }else {
    
    plotting<-"Y"
    
  }


if(plotting=="Y"){
  
  print("Plotting the data...")
  
  #list of genes to annotate
  sign_genes <- setNames(data.frame(unique(sign_pval_DF$gene)), "gene")
  
  
  #######################################################
  #    Annotate genes using BiomaRt if outside the RE   #
  #######################################################
  # Load libraries ----
  if (! requireNamespace("biomaRt", quietly = TRUE)) {
    install.packages("biomaRt")
    loadNamespace("biomaRt")
  }
  
  library("biomaRt")
  
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  annotation<-getBM(attributes=c('hgnc_symbol','chromosome_name', 'start_position', 'end_position','uniprot_gn_symbol','uniprotswissprot','description','strand'),
                    filters=c('hgnc_symbol'),
                    values=list(sign_genes$gene),
                    mart=ensembl)

  annotation <- annotation[!grepl("CHR_.*|KI.*|GL.*", annotation$chromosome_name), ]
  
  
  #Order and remove duplicates
  annotation <- annotation %>% arrange(desc(uniprotswissprot), desc(chromosome_name)) %>% distinct(hgnc_symbol, .keep_all = TRUE)
  
  
  #Clean description
  
  annotation$description <- sub("\\[.*","",annotation$description)

  #Annotate significant genes
  sign_genes = merge(x=sign_genes, y=annotation, by.x=c("gene"), by.y=c("hgnc_symbol"),all.x = TRUE)
  
  
  print(paste0("There are ",nrow(sign_pval_DF)," significant signals with p.adjust.fdr <= ",p.adjust.fdr_threshold, " in ",nrow(sign_genes)," genes"))
  
  #Check missing annotation in significant genes
  missing_annotation<-sign_genes %>% filter(is.na(sign_genes$uniprotswissprot))
  
  if(nrow(missing_annotation)>0){
    
    print(paste0(nrow(missing_annotation), " significant genes out of ",nrow(sign_genes)," have no BioMart annotation: ", paste0(paste0("'", missing_annotation$gene, "'"), collapse = ", ")))
    
  }else{
    
    print("no BioMart annotation missing for significant genes")
    
  }

  #update pval FDR dataframe
  pval_DF_annotated = merge(x=pval_DF, y=sign_genes[,c("gene",colnames(sign_genes)[!colnames(sign_genes) %in% colnames(pval_DF)])], by=c("gene"),all.x = TRUE)
  
  #Check
  sum(is.na(pval_DF_annotated$chromosome_name))
  sum(is.na(pval_DF_annotated$start_position))
  sum(is.na(pval_DF_annotated$uniprotswissprot))

  #Filter significant FDR signals
  sign_pval_DF<-pval_DF_annotated %>% filter(p.adjust.fdr<=p.adjust.fdr_threshold)
  
  
  
  #Extract list of unique disease level per significant signals
  sign_analysis_label <- setNames(as.data.frame(table(sign_pval_DF$analysis.label)), c("analysis.label", "Freq"))
  
  
    setwd(path_to_analysis)
    
    #extract pvalues for disease
    pval_DF_analysis_label<-pval_DF_annotated %>% dplyr::filter(analysis.label==disease)
    
    #Show analysis disease label
    sign_analysis_label_read <- as.data.frame(table(pval_DF_analysis_label$analysis))
    sign_analysis_label_read<-as.character(sign_analysis_label_read$Var1)
    print(sign_analysis_label_read)
    
    #create analysis disease folder
    path<-paste0("plots/",sign_analysis_label_read)
    dir.create(file.path(getwd(), path), recursive = TRUE)
    setwd(file.path(getwd(), path))
    getwd()
    
    #Significant pvalues above threshold
    pval_DF_analysis_label_sign<-pval_DF_analysis_label %>% dplyr::filter(p.adjust.fdr<=p.adjust.fdr_threshold)
    
    print(paste0(nrow(pval_DF_analysis_label_sign), " significant gene-disease associations found in association with ",sign_analysis_label_read))
    
    # Load additional library for plotting ----
    if (! requireNamespace("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
      loadNamespace("ggplot2")
    }
    if (! requireNamespace("ggrepel", quietly = TRUE)) {
      install.packages("ggrepel")
      loadNamespace("ggrepel")
    }

    library(ggplot2)
    library(ggrepel)
    
    ###############################
    #       Volcano plot        #
    ###############################
    
    pval_DF_analysis_label$p.adjust.fdr<-as.numeric(pval_DF_analysis_label$p.adjust.fdr)
      pval_DF_analysis_label$log_fdr<-(-log10(pval_DF_analysis_label$p.adjust.fdr))
    if (any(is.infinite(pval_DF_analysis_label$log_fdr))) {
    if (all(is.infinite(pval_DF_analysis_label$log_fdr))) {
      # Set all log_or values to 5 if all are infinite
      pval_DF_analysis_label$log_fdr <- 8
      print("All -log10(p.adjust.fdr) go to infinite: value set to default")
      
    }else{
      
      pval_DF_analysis_label$log_fdr<-ifelse(is.infinite(pval_DF_analysis_label$log_fdr),max(-log10(pval_DF_analysis_label$p.adjust.fdr[is.finite(-log10(pval_DF_analysis_label$p.adjust.fdr))]))+1,pval_DF_analysis_label$log_fdr)
      print("-log10(p.adjust.fdr) that to go to infinite are approximated to the maximum value")
      
    }
    }
    
    pval_DF_analysis_label$or<-as.numeric(pval_DF_analysis_label$or)
      pval_DF_analysis_label$log_or<-(log2(pval_DF_analysis_label$or))
      if (any(is.infinite(pval_DF_analysis_label$log_or))) {
      if (all(is.infinite(pval_DF_analysis_label$log_or))) {
        # Set all log_or values to 5 if all are infinite
        pval_DF_analysis_label$log_or <- 5
        print("All log2(OR) go to infinite: value set to default")
        
      } else {
        
        pval_DF_analysis_label$log_or<-ifelse(is.infinite(pval_DF_analysis_label$log_or),max(log2(pval_DF_analysis_label$or[is.finite(log2(pval_DF_analysis_label$or))]))+1,pval_DF_analysis_label$log_or)
        print("-log2(OR) that to go to infinite are approximated to the maximum value")
      }
      }
    


      #Generate volcano plot
      
      sign_to_plot<-pval_DF_analysis_label %>% dplyr::filter(p.adjust.fdr<= p.adjust.fdr_threshold & or >= 3)
      
      log_fdr_intercept<-(-log10(p.adjust.fdr_threshold))
      log_or_intercept<-log2(3)
      
      volcano <- ggplot(pval_DF_analysis_label, aes(x = log_or, y = log_fdr, color = test)) +
        geom_point() +
        theme_classic() +  # Use classic theme for black axis lines
        
        # Dashed intercept lines
        geom_hline(yintercept = log_fdr_intercept, col = "maroon", linetype = "dashed") +
        geom_vline(xintercept = log_or_intercept, col = "maroon", linetype = "dashed") +
        
        # Set axis limits to start from 0, add a small extension to avoid cutting off points, and set custom breaks
        scale_x_continuous(limits = c(0, max(pval_DF_analysis_label$log_or) * 1.05), expand = c(0, 0), breaks = scales::pretty_breaks()) +
        scale_y_continuous(limits = c(0, max(pval_DF_analysis_label$log_fdr) * 1.05), expand = c(0, 0), breaks = scales::pretty_breaks()) +
        
        # Add labels and axis titles
        labs(
          title = sign_analysis_label_read,
          x = expression(log[2](OR)),
          y = expression(-log[10](p.adjust.fdr))
        ) +
        
        # Add gene labels to specific points
        geom_text_repel(data = sign_to_plot, size = 4, aes(label = gene), show.legend = FALSE, colour = 'gray20') +
        
        # Customize theme to reduce space between axis and label and add minor grid lines
        theme(
          text = element_text(family = "sans"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(margin = margin(t = 5)),  # Reduce space for x-axis label
          axis.title.y = element_text(margin = margin(r = 5)),  # Reduce space for y-axis label
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Add smaller grid lines
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          plot.title = element_text(
            size = 16,
            face = "bold",
            family = "sans",
            color = "gray60",
            hjust = 0.5,
            lineheight = 1.2
          )
        )
      
      # Save the plot
      ggsave(volcano, file = paste0(disease, "_volcanoplot.png"), width = 2560/227, height = 1600/227, dpi = 227, bg = "#ffffff")
      
      # Print save confirmation
      print(paste0(disease, "_volcanoplot.png saved!"))
      

      
    ###############Generate addiyional tables,lolliplot and hpo plot
  
    #read caco definition TSV file
  caco_def<-read_tsv(paste0(path_to_analysis,"/data/",disease,".tsv"))

  print(paste0("Reading ",disease,".tsv file"))
  
  #Create loop per each significant gene
  for(sign_gene in unique(pval_DF_analysis_label_sign$gene)){
    
    print(sign_gene)
    
    #Get significant pvalues per that gene
    pval_sign_gene<-pval_DF_analysis_label_sign %>% dplyr::filter(gene==sign_gene)
    
    #Create loop per each significant test
    for(sign_test in pval_sign_gene$test){
      
      sign_gene_test<-paste0(sign_gene,"_",sign_test)
      
      print(sign_gene_test)
      
      #Create folder for each significan gene-test
      dir.create(path = paste0(path_to_analysis,"/plots","/",sign_analysis_label_read,"/",sign_gene_test), recursive = TRUE)
      path<-paste0(path_to_analysis,"/plots","/",sign_analysis_label_read,"/",sign_gene_test)
      setwd(path)
      getwd()
      
      #Extract row for significant gene-test
      pval_sign_gene_test<-pval_sign_gene %>% dplyr::filter(gene==sign_gene) %>% dplyr::filter(test==sign_test)
      
      #Extract all variants in that gene
      gene_all<-exomiser  %>% dplyr::filter(gene==sign_gene)
      
      #Extract caco or caco.denovo info
      if(sign_test=="denovo"){
        
        gene_all$caco<-caco_def$caco.denovo[match(gene_all$sample.id,caco_def$sample.id)]
        
        
      }else{
        
        gene_all$caco<-caco_def$caco[match(gene_all$sample.id,caco_def$sample.id)]
        
      }
      
      
      table(gene_all$caco, useNA = "always")
      
      
      if (file.exists(paste0(path_to_analysis, "/output/probandDataMatrix/probandDataMatrix_", disease, ".tsv"))) {
        
        #####################
        contrib_variants <- read_tsv(paste0(path_to_analysis, "/output/probandDataMatrix/probandDataMatrix_", disease, ".tsv"))
        
        print(paste0("Reading probandDataMatrix_", disease, ".tsv file"))
        
        contrib_variants <- contrib_variants %>%
          dplyr::filter((drop_1 == 1 & drop_2 == 0) | (drop_1 == 0 & drop_2 == 1) | (drop_1 == 0 & drop_2 == 0) | (drop_1 == 0 & is.na(drop_2))) %>%
          dplyr::filter(gene == sign_gene)
        
        # Filter variants based on the test that has been run
        if (sign_test == "LoF") {
          contrib_variants <- contrib_variants %>% dplyr::filter(max.vartype == "LoF")
        } else if (sign_test == "zero80") {
          contrib_variants <- contrib_variants %>% dplyr::filter(max.varscore >= 0.80)
        } else if (sign_test == "denovo") {
          contrib_variants <- contrib_variants %>% dplyr::filter(max.denovo == '1' & fam.structure >= 3)
        } else {
          contrib_variants <- contrib_variants %>% dplyr::filter(max.ccr == 1 & max.varscore >= 0.8)
        }
        
        # Extract caco or caco.denovo info
        if (sign_test == "denovo") {
          contrib_variants$caco <- caco_def$caco.denovo[match(contrib_variants$sample.id, caco_def$sample.id)]
        } else {
          contrib_variants$caco <- caco_def$caco[match(contrib_variants$sample.id, caco_def$sample.id)]
        }
        
        table(contrib_variants$caco, useNA = "always")
        
        cases_contrib <- contrib_variants %>% dplyr::filter(caco == 1)
        controls_contrib <- contrib_variants %>% dplyr::filter(caco == 0)
        
        if (length(unique(cases_contrib$sample.id)) == pval_sign_gene_test$d) {
          print(paste0("Match: ", length(unique(cases_contrib$sample.id)), " variants found == d (", pval_sign_gene_test$d, ")"))
        } else {
          print(paste0("Mismatch: ", length(unique(cases_contrib$sample.id)), " variants found != d (", pval_sign_gene_test$d, ")"))
          break
        }
        
        if (length(unique(controls_contrib$sample.id)) == pval_sign_gene_test$c) {
          print(paste0("Match: ", length(unique(controls_contrib$sample.id)), " variants found == c (", pval_sign_gene_test$c, ")"))
        } else {
          print(paste0("Mismatch: ", length(unique(controls_contrib$sample.id)), " variants found != c (", pval_sign_gene_test$c, ")"))
          break
        }
        
        contrib_df <- contrib_variants %>%
          dplyr::select(sample.id, variant.lift_1, variant.lift_2) %>%
          reshape2::melt(id.vars = "sample.id") %>%
          dplyr::filter(!is.na(value))
        
      }
      
      #Add readable caco definition
      gene_all$caco_definition<-ifelse(gene_all$caco==1,"case",NA)
      gene_all$caco_definition<-ifelse((is.na(gene_all$caco_definition)==TRUE) & gene_all$caco==0,"control",gene_all$caco_definition)
      gene_all$caco_definition<-ifelse((is.na(gene_all$caco_definition)==TRUE),"excluded",gene_all$caco_definition)
      
      table(gene_all$caco_definition, useNA = "always")
      
      #Filter variants based on the test that has been run
      if (sign_test=="LoF"){
        gene_test<- gene_all  %>% dplyr::filter(vartype=="LoF")
      }else if(sign_test=="zero80"){
        gene_test<- gene_all  %>% dplyr::filter(variant.score>=0.80)
      }else if(sign_test=="denovo"){
        gene_test<- gene_all  %>% dplyr::filter(de.novo=='Y' & fam.structure >=3)
      }else{
        gene_test<- gene_all  %>% dplyr::filter(ccr.flag==1 & variant.score>=0.8)
      }
      
      if (exists("contrib_df")) {
        
      gene_test$check <- ifelse(gene_test$sample.id %in% contrib_df$sample.id & gene_test$variant.lift %in% contrib_df$value, "Contrib", "Not contrib")
      
      gene_test$caco_definition <- ifelse(gene_test$caco_definition == "case" & gene_test$check == "Not contrib", "case_drop_varfilter", gene_test$caco_definition)
      gene_test$caco_definition <- ifelse(gene_test$caco_definition == "control" & gene_test$check == "Not contrib", "control_drop_varfilter", gene_test$caco_definition)
      gene_test$check<-NULL
      
      }
      
      table(gene_test$caco_definition, useNA = "always")
      
      # Extract number of cases contributing to the test
      cases <- gene_test %>% dplyr::filter(caco_definition == "case")
      controls <- gene_test %>% dplyr::filter(caco_definition == "control")
      
      if(length(unique(cases$sample.id))==pval_sign_gene_test$d){
        
        print(paste0("Match: ",length(unique(cases$sample.id))," variants found == d (",pval_sign_gene_test$d,")"))
        
      }else{
        
        print(paste0("Mismatch: ",length(unique(cases$sample.id))," variants found != d (",pval_sign_gene_test$d,")"))
        
      }
      
      
      if (length(unique(controls$sample.id)) == pval_sign_gene_test$c) {
        print(paste0("Match: ", length(unique(controls$sample.id)), " variants found == c (", pval_sign_gene_test$c, ")"))
      } else {
        print(paste0("Mismatch: ", length(unique(controls$sample.id)), " variants found != c (", pval_sign_gene_test$c, ")"))
        break
      }
      
      # Extract not contributing variants based on the test that has been run
      if (sign_test == "LoF") {
        not_contributing <- gene_all %>% dplyr::filter(vartype != "LoF")
      } else if (sign_test == "zero80") {
        not_contributing <- gene_all %>% dplyr::filter(variant.score < 0.80)
      } else if (sign_test == "denovo") {
        not_contributing <- gene_all %>% dplyr::filter(de.novo != 'Y' | fam.structure < 3)
      } else {
        not_contributing <- gene_all %>% dplyr::filter(ccr.flag != 1 | variant.score < 0.8)
      }
      
      not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "case", "not_contributing_case", not_contributing$caco_definition)
      not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "control", "not_contributing_control", not_contributing$caco_definition)
      not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "excluded", "not_contributing_excluded", not_contributing$caco_definition)

      #Merge contributing and not contributing variants
      gene_test<-rbind(gene_test,not_contributing)
      
      #Add additional info
      gene_test <-merge(x=gene_test, y=pval_sign_gene_test[,c("gene",colnames(pval_sign_gene_test)[!colnames(pval_sign_gene_test) %in% colnames(gene_test)])], by=c("gene"),all.x=TRUE)
      
      #Order columns
      gene_test<-gene_test %>% dplyr::select(caco_definition,caco,sample.id,hpo.terms,hpo.ids,gene,fam.structure,variant:ccr.flag,vartype,everything()) %>% dplyr::arrange(caco_definition,variant)
      
      #Eventually do some data cleaning
      gene_test$hgvs<-sub("\\|$","",gene_test$hgvs)
      gene_test$exomiser.result.count<-sub("\\|$","",gene_test$exomiser.result.count)
      gene_test$genotypes<-sub("\\|$","",gene_test$genotypes)

      
      # Order table by case/control definition
      caco_def_order <- c("case", "case_drop_varfilter", "excluded", "control", "control_drop_varfilter", 
                          "not_contributing_case", "not_contributing_excluded", "not_contributing_control")
      
      gene_test <- gene_test[order(factor(gene_test$caco_definition, levels = unique(caco_def_order))), ]
      
      #Write variants table
      write.table(gene_test,paste0(sign_gene_test,"_",disease,"_variants.tsv"), sep="\t", row.names = FALSE)
      
      print(paste0(sign_gene_test,"_",disease,"_variants.tsv written!"))
      
      
      ###############################
      #           Lolliplot         #
      ###############################
      
      
      #Plot only cases with contruting variants
      gene_test<-gene_test %>% dplyr::filter(caco_definition=="case")
      
      #Impose limit of 100 variants to be plotted, plot wouldn't be readable
      if(nrow(gene_test)>=100){
        
        print(paste0("Error: ",nrow(gene_test), " variants found in cases and exluded probands, visualisation limit set to 100"))
        
        next
      }
      
      #re-annotate genotype
      gene_test$genotype <- ifelse(gene_test$exomisermoi=="AD",'het',NA)
      
      pattern<-".*\\|.*"
      gene_test$genotype <- ifelse(is.na(gene_test$genotype)==TRUE & grepl(pattern,gene_test$functional.class),'comp_het',gene_test$genotype)
      
      gene_test$genotype <- ifelse(is.na(gene_test$genotype)==TRUE & gene_test$exomisermoi=="AR",'hom',gene_test$genotype)
      
      gene_test$genotype <- ifelse(is.na(gene_test$genotype)==TRUE & gene_test$exomisermoi=="XD",'x_het',gene_test$genotype)
      
      gene_test$genotype <- ifelse(is.na(gene_test$genotype)==TRUE & gene_test$exomisermoi=="XR",'x_hom',gene_test$genotype)
      
      sum(is.na(gene_test$genotype))
      table(gene_test$genotype,useNA = "always")
      
      #Order by mode of inheritance
      gene_test$exomisermoi <- factor(gene_test$exomisermoi, levels = c("AR", "AD","XR","XD"))
      
      gene_test <- gene_test[order(gene_test$exomisermoi),]
      
      #split comp het variants
      gene_test <- separate_rows(gene_test, sample.id, sep = "\\|") %>%
        separate_rows(c('variant','functional.class','hgvs','exomiser.result.count'), sep = "\\|")
      
      #split HGVS info
      gene_test<-separate(data = gene_test, col = hgvs, into = c("hgvs_gene_id", "hgvs_transcript","hgvs_c_change","hgvs_p_change"), sep = ":")
      
      #remove duplicates
      gene_test<-gene_test[!duplicated(gene_test[,c('sample.id','gene','variant')]),]
      
      #info abou the gene
      gene_info<-pval_sign_gene %>% dplyr::slice(1) %>% dplyr::select(gene,chromosome_name,start_position,end_position,description,strand,uniprotswissprot)
      
      #Uniprot code for the gene
      uniprot<-gene_info %>% dplyr::select(gene,uniprotswissprot)
      
      if(is.na(uniprot$uniprotswissprot)){
        
        next
        print(paste0("no uniprot id for ", sign_gene, " gene!"))
        
      }
      
      #Data wrangling
      gene_test$hgvs_p_change <- sub("p.","",gene_test$hgvs_p_change)
      gene_test$hgvs_p_change <- sub("[(]","",gene_test$hgvs_p_change)
      gene_test$hgvs_p_change <- sub("[)]","",gene_test$hgvs_p_change)

      # Load additional library for lolliplot ----
      if (! requireNamespace("httr", quietly = TRUE)) {
        install.packages("httr")
        loadNamespace("httr")
      }
      if (! requireNamespace("drawProteins", quietly = TRUE)) {
        install.packages("drawProteins")
        loadNamespace("drawProteins")
      }

      library(httr)
      library(drawProteins)

      #Retrieve protein info
      protein_id<-uniprot[uniprot$gene==sign_gene,'uniprotswissprot']
      protein_id_url<-gsub(" ","%2C",protein_id)
      baseurl<-"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
      url<-paste0(baseurl,protein_id_url)
      prots_feat<-GET(url,accept_json())
      prots_feat_red<-httr::content(prots_feat)
      features_total_plot<-NULL
      
      for(k in 1:length(prots_feat_red)){ 
        features_temp<-drawProteins::extract_feat_acc(prots_feat_red[[k]])#the extract_feat_acc() function takes features into a data.frame
        features_temp$order<-k # this order is needed for plotting later
        features_total_plot<-rbind(features_total_plot,features_temp)
      }
      
      plot_start<-0#starts 0
      plot_end<-max(features_total_plot$end,na.rm = TRUE)
      
      #Break if variants are not all on b38 and not all annotated on same transcript
      if(all(gene_test$assembly=="b38")){
        
        
        transcripts<-as.data.frame(table(gene_test$hgvs_transcript))
        
        if(nrow(transcripts)>1){
          
          print(paste0("Variants are NOT all annotated on the same transcipt: ", paste0(paste0("'", transcripts$Var1, "'"), collapse = ", ")))
          
          break
          
        } else{
          
          print(paste0("Variants are all annotated on the same transcipt: ", paste0(paste0("'", transcripts$Var1, "'"), collapse = ", ")))
          
          
        }
        
        
      
        transcripts<-as.data.frame(table(gsub("\\..*","",transcripts$Var1)))
        
        select_transcript<-as.character(transcripts$Var1)
        
      }else{
        
        print("Not all variants are annotated on b38 MANE transcript")
        
        check_assembly<- gene_test %>% dplyr::filter(!assembly=="b38")
        
        break
        
      }
      
      
      gene_test$select_transcript<-as.character(select_transcript)
      
      
      #Re-annotate variants with no HGVS annotated protein change by identifying closer aminoacid
      gene_test$select_p_change<-NA
      
      gene_test$select_p_change<-ifelse((grepl(paste0("^",select_transcript),gene_test$hgvs_transcript))&((gene_test$hgvs_p_change=='?'|gene_test$hgvs_p_change=='=')==FALSE),gene_test$hgvs_p_change,gene_test$select_p_change)
      
      fix_transcript<-gene_test %>% dplyr::filter(is.na(select_p_change))
      fix_transcript<-fix_transcript[!duplicated(fix_transcript[,c('gene','hgvs_c_change','assembly')]),]
      
      print(paste0(nrow(fix_transcript)," variants to fix!"))
      
      fix_transcript$fixed<-"Y"
      
      if((nrow(fix_transcript)>0)==TRUE){
        
        # Load additional library for re-annotation ----
        if (! requireNamespace("ensembldb", quietly = TRUE)) {
          install.packages("ensembldb")
          loadNamespace("ensembldb")
        }
        if (! requireNamespace("AnnotationHub", quietly = TRUE)) {
          install.packages("AnnotationHub")
          loadNamespace("AnnotationHub")
        }
        library(ensembldb)
        library(AnnotationHub)
        
        AH <- AnnotationHub()
        
        query(AH, "Hsapiens.v109")
        
        edb <- AH[["AH109606"]]      
        
        edb <- filter(edb, filter = GeneNameFilter(sign_gene))
        
        
        for (k in 1:length(fix_transcript$gene)){
          tryCatch({
            w<-fix_transcript[k,]
            print(w$hgvs_c_change)
            
            
            test<-separate(data = w, col = variant, into = c("chromosome", "var","wt","mut"), sep = ":",remove = FALSE)
            test$var<-as.numeric(test$var)
            var<-test$var
            print(var)
            
            suppressWarnings({
              
              edb <- filter(edb, filter = ~ seq_name == test$chromosome)
              gnm_pos <- GRanges(test$chromosome, IRanges(var, width = 1))
              trc_pos <- genomeToTranscript(gnm_pos, edb)
              t_data<-data.frame(trc_pos)
              
              if(any(colnames(t_data)=='names')){
                
                t_data<-t_data %>% dplyr::filter(names==gsub("\\..*","",select_transcript))
                
                if(nrow(t_data)>0){
                  
                  rng_tx <- IRanges(start = t_data$start , width = t_data$width,
                                    names = t_data$names)
                  rng_prt <- transcriptToProtein(rng_tx, edb)
                  p_data<-data.frame(rng_prt)
                  fix_transcript$select_p_change[k]<-p_data$start
                  
                }
              }
              
            })
            
            if(grepl("^c.\\-",fix_transcript$hgvs_c_change[k])==TRUE){
              
              
              fix_transcript$select_p_change[k]<-'1'
              
            }
            
            
            counter<-0
            while(((fix_transcript$select_p_change[k]==-1)|(is.na(fix_transcript$select_p_change[k])))&(counter<3000)){
              counter <- sum(counter, 1)
              
              if(grepl("c.\\*",fix_transcript$hgvs_c_change[k])==TRUE){
                fix_transcript$select_p_change[k]<-"end"
                break
              }else if(grepl("^c.\\-",fix_transcript$hgvs_c_change[k])==TRUE){
                
                var<-ifelse(((grepl("^c.\\-",fix_transcript$hgvs_c_change[k])==TRUE)&(gene_info$strand==1)),(var+1),(var-1))
                print(var)
                suppressWarnings({
                  edb <- filter(edb, filter = ~ seq_name == test$chromosome)
                  gnm_pos <- GRanges(test$chromosome, IRanges(var, width = 1))
                  trc_pos <- genomeToTranscript(gnm_pos, edb)
                  t_data<-data.frame(trc_pos)
                  if(any(colnames(t_data)=='names')){
                    
                    t_data<-t_data %>% dplyr::filter(names==gsub("\\..*","",select_transcript))
                    
                    if(nrow(t_data)>0){
                      
                      rng_tx <- IRanges(start = t_data$start , width = t_data$width,
                                        names = t_data$names)
                      rng_prt <- transcriptToProtein(rng_tx, edb)
                      p_data<-data.frame(rng_prt)
                      fix_transcript$select_p_change[k]<-p_data$start
                      
                    }
                  }
                })
                
              }else{
                word<-c("\\+")
                pattern <- paste0('.*', word, '.*', collapse = '|')
                var<-ifelse(((grepl(pattern,fix_transcript$hgvs_c_change[k])==TRUE)&(gene_info$strand==-1))|((grepl(pattern,fix_transcript$hgvs_c_change[k])==FALSE)&(gene_info$strand==1)&(grepl("_",fix_transcript$hgvs_c_change[k])==FALSE)),(var+1),(var-1))
                print(var)
                suppressWarnings({
                  edb <- filter(edb, filter = ~ seq_name == test$chromosome)
                  gnm_pos <- GRanges(test$chromosome, IRanges(var, width = 1))
                  trc_pos <- genomeToTranscript(gnm_pos, edb)
                  t_data<-data.frame(trc_pos)
                  if(any(colnames(t_data)=='names')){
                    
                    t_data<-t_data %>% dplyr::filter(names==gsub("\\..*","",select_transcript))
                    
                    if(nrow(t_data)>0){
                      
                      rng_tx <- IRanges(start = t_data$start , width = t_data$width,
                                        names = t_data$names)
                      rng_prt <- transcriptToProtein(rng_tx, edb)
                      p_data<-data.frame(rng_prt)
                      fix_transcript$select_p_change[k]<-p_data$start
                      
                    }
                  }
                })
                
              }
              
            }
            
            
            var<-test$var
            
            counter<-0
            while(((fix_transcript$select_p_change[k]==-1)|(is.na(fix_transcript$select_p_change[k])))&(counter<3000)){
              
              counter <- sum(counter, 1)
              
              if(grepl("c.\\*",fix_transcript$hgvs_c_change[k])==TRUE){
                fix_transcript$select_p_change[k]<-"end"
                break
              }else{
                word<-c("\\+")
                pattern <- paste0('.*', word, '.*', collapse = '|')
                var<-ifelse(((grepl(pattern,fix_transcript$hgvs_c_change[k])==TRUE)&(gene_info$strand==-1))|((grepl(pattern,fix_transcript$hgvs_c_change[k])==FALSE)&(gene_info$strand==1)&(grepl("_",fix_transcript$hgvs_c_change[k])==FALSE)),(var-1),(var+1))
                print(var)
                suppressWarnings({
                  edb <- filter(edb, filter = ~ seq_name == test$chromosome)
                  gnm_pos <- GRanges(test$chromosome, IRanges(var, width = 1))
                  trc_pos <- genomeToTranscript(gnm_pos, edb)
                  t_data<-data.frame(trc_pos)
                  if(any(colnames(t_data)=='names')){
                    
                    t_data<-t_data %>% dplyr::filter(names==gsub("\\..*","",select_transcript))
                    
                    if(nrow(t_data)>0){
                      
                      rng_tx <- IRanges(start = t_data$start , width = t_data$width,
                                        names = t_data$names)
                      rng_prt <- transcriptToProtein(rng_tx, edb)
                      p_data<-data.frame(rng_prt)
                      fix_transcript$select_p_change[k]<-p_data$start
                      
                    }
                  }
                })
                
                
              }
              
            }
            
            var<-test$var
            
            while((fix_transcript$select_p_change[k]==-1)|(is.na(fix_transcript$select_p_change[k]))){
              
              if(grepl("c.\\*",fix_transcript$hgvs_c_change[k])==TRUE){
                fix_transcript$select_p_change[k]<-"end"
                break
              }else{
                word<-c("\\+")
                pattern <- paste0('.*', word, '.*', collapse = '|')
                var<-ifelse(((grepl(pattern,fix_transcript$hgvs_c_change[k])==TRUE)&(gene_info$strand==-1))|((grepl(pattern,fix_transcript$hgvs_c_change[k])==FALSE)&(gene_info$strand==1)&(grepl("_",fix_transcript$hgvs_c_change[k])==FALSE)),(var+100),(var-100))
                print(var)
                suppressWarnings({
                  edb <- filter(edb, filter = ~ seq_name == test$chromosome)
                  gnm_pos <- GRanges(test$chromosome, IRanges(var, width = 1))
                  trc_pos <- genomeToTranscript(gnm_pos, edb)
                  t_data<-data.frame(trc_pos)
                  if(any(colnames(t_data)=='names')){
                    
                    t_data<-t_data %>% dplyr::filter(names==gsub("\\..*","",select_transcript))
                    
                    if(nrow(t_data)>0){
                      
                      rng_tx <- IRanges(start = t_data$start , width = t_data$width,
                                        names = t_data$names)
                      rng_prt <- transcriptToProtein(rng_tx, edb)
                      p_data<-data.frame(rng_prt)
                      fix_transcript$select_p_change[k]<-p_data$start
                      
                    }
                  }
                })
                
              }
              
            }
            skip_to_next <- FALSE
          }, error=function(e) {cat('Error:no uniprot id!');skip_to_next <<- TRUE})
          if(skip_to_next) { next }
        }
        
      }
      
      
      gene_test$select_p_change<-ifelse(is.na(gene_test$select_p_change),fix_transcript$select_p_change[match(gene_test$hgvs_c_change,fix_transcript$hgvs_c_change)],gene_test$select_p_change)
      
      gene_test$select_p_change<-ifelse(gene_test$select_p_change=="end",as.character(plot_end),gene_test$select_p_change)
      
      gene_test$fixed<-fix_transcript$fixed[match(gene_test$select_p_change,fix_transcript$select_p_change)]

      gene_test$hgvs_p_change<-ifelse((gene_test$hgvs_p_change=='?'|gene_test$hgvs_p_change=='='),"",gene_test$hgvs_p_change)
      
      #Show re-annotated variants between parenthesis ( )
      gene_test$fixed_p_change<-ifelse(((gene_test$hgvs_p_change==gene_test$select_p_change)==FALSE)&(gene_test$fixed=="Y"),paste0(gene_test$hgvs_p_change,"\n(",gene_test$select_p_change,")"),gene_test$select_p_change)
      
      gene_test$fixed_p_change<-ifelse(is.na(gene_test$fixed_p_change),gene_test$select_p_change,gene_test$fixed_p_change)
      
      #Rename two columns
      gene_test <- gene_test %>% dplyr::rename(gene.symbol = hgvs_gene_id, protein.change = select_p_change)

      
      #Create amino acid number variable
      var.aanum<-unlist(lapply(regmatches(gene_test$protein.change,gregexpr(pattern="*(\\d+)",gene_test$protein.change)),function(x) x[[1]]))
      gene_test<-cbind(gene_test,as.numeric(as.character(var.aanum)))#to remove factor
      colnames(gene_test)[ncol(gene_test)]<-"var.aanum"
      
      #Define end of plot
      plot_end<-ifelse((max(gene_test$var.aanum)>plot_end)==TRUE,max(gene_test$var.aanum),plot_end)
      
      gene_test$protein.change<-as.character(gene_test$protein.change)
      
      gene_test$protein.change<-ifelse(gene_test$protein.change=="end",as.character(plot_end),gene_test$protein.change)
      
      #Additional annotation of the variants
      var<-gene_test
      var_allvariant<-var %>% mutate_all(as.character) #%>% select(gene,sample.id,variant,protein.change,functional.class,genotype,HL.def.hpo,syndromic,other.def.hpo,other.def.ancestors)
      var<-var %>% dplyr::select(caco_definition,gene,sample.id,protein.change,functional.class,genotype,fixed_p_change,hpo.terms,hpo.ids)
      var<-unique(var[!duplicated(var),])#remove duplicates
      var<-var[order(var$gene,gsub("([A-Z]+)([0-9]+)","\\2",var$protein.change)),]#order
      var.freq<-var %>%
        dplyr::group_by(fixed_p_change,functional.class,genotype) %>%
        dplyr::mutate(freq = n()) %>% 
        dplyr::slice(1)
      
      var.freq<-var.freq %>% ungroup 
      
      #Create amino acid number variable
      var.aanum<-unlist(lapply(regmatches(var.freq$protein.change,gregexpr(pattern="*(\\d+)",var.freq$protein.change)),function(x) x[[1]]))
      var.plot.data<-cbind(var.freq,as.numeric(as.character(var.aanum)))#hard way to remove factor
      colnames(var.plot.data)[ncol(var.plot.data)]<-"var.aanum"

      #Data wrangling before plotting
      var.plot.table<-var.plot.data[var.plot.data$gene==sign_gene,]
      var.plot.table<-var.plot.table[order(var.plot.table$var.aanum),]
      
      var.plot<-var.plot.data[var.plot.data$gene==sign_gene,]
      var.plot<-var.plot[order(var.plot$var.aanum),]
      var.plot$protein.change<-as.character(var.plot$protein.change)
      var.plot<-var.plot %>% mutate_all(as.character)
      var.plot$freq<-as.numeric(var.plot$freq)
      var.plot$var.aanum<-as.numeric(var.plot$var.aanum)
      var.plot$genotype <- factor(var.plot$genotype, levels = c("het", "x_het","hom","x_hom","comp_het"))
      var.plot <- var.plot[order(var.plot$genotype),]
      
      #Create shapeflag for different genotype
      var.plot$shapeflag <- ifelse((var.plot$genotype=="het"|var.plot$genotype=="x_het"), "16", NA)
      var.plot$shapeflag <- ifelse(is.na(var.plot$shapeflag)==TRUE & (var.plot$genotype=="hom"|var.plot$genotype=="x_hom"), "15", var.plot$shapeflag)
      var.plot$shapeflag <- ifelse(is.na(var.plot$shapeflag)==TRUE & var.plot$genotype=="comp_het", "18", var.plot$shapeflag)

      shapeflags<-as.vector(unique(var.plot$shapeflag))
      labelflags<-as.vector(unique(var.plot$genotype))
      
      #Create lolliplot
      p<-ggplot(var.plot,aes(var.aanum,freq,color=functional.class,label=fixed_p_change,group=shapeflag,shape=shapeflag))+
        geom_segment(aes(x=var.aanum,y=0,xend=var.aanum,yend=freq),color=ifelse((var.plot$caco_definition=="excluded"),"gold","grey50"),linewidth=ifelse(var.plot$freq>=1,0.7,0.7))+
        geom_point(aes(shape=shapeflag),size=5.5)+
        scale_shape_manual(values=c('16'=16,'15'=15,'18'=18),labels=c('16'="het",'15'="hom",'18'="comp_het"),name="Genotypes")+
        geom_text_repel(nudge_y=0.2,color=ifelse(var.plot$freq>=1,"black","NA"),size=ifelse(var.plot$freq>=1,4,4),fontface="bold",max.overlaps = Inf)
      
      #Extract and include additional protein annotation
      features_total_plot<-features_total_plot %>% ungroup()
      
      features_total_plot<-merge(features_total_plot, data.frame(table(type = features_total_plot$type)), by = "type")
      features_total_plot<-features_total_plot %>% dplyr::filter(Freq<=8) %>% dplyr::filter(!description=="Disordered")
      
      if("CHAIN"%in%(unique(features_total_plot$type))) p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="CHAIN",],mapping=aes(xmin=plot_start,xmax=plot_end,ymin=-0.45,ymax=0),colour="grey85",fill="grey85",inherit.aes=F)
      if("DNA_BIND"%in%(unique(features_total_plot$type))) {
        features_total_plot[features_total_plot$type=="DNA_BIND","description"]<-features_total_plot[features_total_plot$type=="DNA_BIND","type"]
        p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="DNA_BIND",],mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),inherit.aes=F)
      }
      if("DOMAIN"%in%(unique(features_total_plot$type))) p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="DOMAIN",],mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),inherit.aes=F)
      if("REGION"%in%(unique(features_total_plot$type))) p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="REGION",],mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),inherit.aes=F)
      if("MOTIF"%in%(unique(features_total_plot$type))) p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="MOTIF",],mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),inherit.aes=F)
      
      if(any(features_total_plot$type=="DOMAIN")==TRUE){
        features_total_plot<-subset(features_total_plot,type%in% c("CHAIN","DNA_BIND","DOMAIN","REGION","MOTIF"))
      }else{
        if("REPEAT"%in%(unique(features_total_plot$type))) p<-p+geom_rect(data=features_total_plot[features_total_plot$type=="REPEAT",],mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),inherit.aes=F)
        features_total_plot<-subset(features_total_plot,type%in% c("CHAIN","DNA_BIND","DOMAIN","REPEAT","REGION","MOTIF"))
      }
      
      #Plot lolliplot with protein annotation
      p<-p+
        scale_x_continuous(breaks=round(seq(plot_start,plot_end,by=100)),name="Amino acid number")+
        theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),axis.line = element_line(colour = "gray70"),plot.title=element_text(face="bold",size=(15),hjust=0),axis.text.x=element_text(angle=90,hjust=1),plot.subtitle=element_text(size=(13),hjust=0),legend.text=element_text(size=9))+
        labs(y="Mutation frequency in our dataset",title=paste0(sign_gene," variants within the ", sign_analysis_label_read, " cohort"),subtitle=paste0(gene_info$description,"\nprotein structure source: Uniprot(",uniprot$uniprotswissprot,")","\ntranscript: ",select_transcript))+
        scale_y_continuous(breaks=seq(0,max(var.plot$freq),1),name="frequency")
      
      
      ggsave(p, file=paste0(getwd(),"/",sign_gene_test,"_lolliplot_",disease,".png"),width = 20, height = 5, units="in")
      print(paste0(sign_gene_test,"_lolliplot_",disease,".png saved!"))
      
      
      #Save lolliplot table
      text_table<-var_allvariant %>% dplyr::filter(gene==sign_gene)
      
      text_table_unique<-text_table[!duplicated(text_table[,c('sample.id')]),]
      text_table_unique$Patient_ID<-paste0("patient_",rownames(text_table_unique))
      text_table$Patient_ID<-text_table_unique$Patient_ID[match(text_table$sample.id,text_table_unique$sample.id)]
      text_table<-text_table %>% dplyr::select(Patient_ID,everything(),-fixed_p_change)
      write.table(text_table,paste0(getwd(),"/",sign_gene_test,"_lolliplot_",disease,".tsv"), quote = F,sep="\t",row.names = FALSE)

      
      ###############################
      #      HPO plots for cases    #
      ###############################
      
      var_allvariant_plot<-var_allvariant %>% dplyr::filter(caco_definition=="case")
      
      var_allvariant_plot<-var_allvariant_plot[!duplicated(var_allvariant_plot[,c('sample.id')]),]
      
      # Load additional library for HPO analysis ----
      if (! requireNamespace("ontologyIndex", quietly = TRUE)) {
        install.packages("ontologyIndex")
        loadNamespace("ontologyIndex")
      }
      if (! requireNamespace("ontologyPlot", quietly = TRUE)) {
        install.packages("ontologyPlot")
        loadNamespace("ontologyPlot")
      }
      
      library(ontologyIndex)
      library(ontologyPlot)
      
      
      hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo")
      
      ## hpo ids and descriptions
      
      hpo.description <- data.frame(hpo.term=as.character(row.names(as.data.frame(hpo[[2]]))),
                                    hpo.description=as.character(as.data.frame(hpo[[2]])[,1]),stringsAsFactors = F)
      
      
      HPO_case_cohort_unique<-data.frame(table(unlist(strsplit(na.omit(unlist(as.character(var_allvariant_plot$hpo.ids))), ","))))
      HPO_case_cohort_unique$hpo.description<-hpo.description$hpo.description[match(HPO_case_cohort_unique$Var1,hpo.description$hpo.term)]
      HPO_case_cohort_unique$Var1<-as.character(HPO_case_cohort_unique$Var1)
      HPO_case_cohort_unique<-HPO_case_cohort_unique %>% dplyr::arrange(desc(Freq)) 
      write.table(HPO_case_cohort_unique,paste0("./",sign_gene_test,"_hpo_table_cases.tsv"), quote = F,sep="\t",row.names = FALSE)
      print(paste0(sign_gene_test,"_hpo_table_cases.tsv written!"))
      
      
      list_hpo_patient<- list()
      for(patient in 1:length(var_allvariant_plot$sample.id)){
        hpo_terms<-var_allvariant_plot %>% dplyr::slice(patient)
        hpo_terms<-data.frame(table(unlist(strsplit(na.omit(unlist(as.character(hpo_terms$hpo.ids))), ","))))
        hpo_terms<-as.vector(hpo_terms$Var1)
        list_hpo_patient[[paste0("patient_",patient)]] <- hpo_terms
      }
      
      list_hpo_patient
      
      #color by frequency 
      
      colour_by_frequency <- function(
    ontology, 
    terms, 
    term_sets, 
    colour_func=colorRampPalette(c("azure2","slategray3"))
      ) {
        ancestors.by.patient <- lapply(term_sets, function(x) get_ancestors(ontology, x))
        
        patients.with.term.count <- sapply(
          terms,
          function(term) sum(
            sapply(ancestors.by.patient, function(ancs) term %in% ancs)
          )
        )
        
        node.colours <- colour_func(1+diff(range(patients.with.term.count)))[
          patients.with.term.count-min(patients.with.term.count)+1
        ]
        names(node.colours) <- terms
        
        node.colours <- node.colours[which(names(node.colours) %in% terms)]
        
        node.colours
      }
      
      if(nrow(var_allvariant_plot)>1 & nrow(var_allvariant_plot)<4){
        
        
        hpo_plot<-onto_plot(hpo, term_sets = list_hpo_patient, label=label_by_term_set,widht = 1,
                            fillcolor = colour_by_frequency,
                            shape = "circle",fontsize = 45,
                            fixedsize = "true")
        
        hpo_plot
        
        #write_dot(hpo_plot,file=paste0("hpo_plot_",i))
        print(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png saved!"))
        png(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png"), width = 2400, height =1200) 
        plot(hpo_plot)
        dev.off()
        
        
        
      } else if(nrow(var_allvariant_plot)>=4 & nrow(var_allvariant_plot)<7){
        
        
        
        
        hpo_plot<-onto_plot(hpo, term_sets = list_hpo_patient, label=label_by_term_set,widht = 1,
                            fillcolor = colour_by_frequency,
                            shape = "circle",fontsize = 70,
                            fixedsize = "true")
        
        hpo_plot
        
        #write_dot(hpo_plot,file=paste0("hpo_plot_",i))
        print(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png saved!"))
        png(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png"), width = 2400, height =1200) 
        plot(hpo_plot)
        dev.off()
        
        
      }else if(nrow(var_allvariant_plot)>=7 & nrow(var_allvariant_plot)<15){
        
        
        
        hpo_plot<-onto_plot(hpo, term_sets = list_hpo_patient, label=label_by_term_set,widht = 1,
                            fillcolor = colour_by_frequency,
                            shape = "circle",fontsize = 100,
                            fixedsize = "true")
        
        hpo_plot
        
        #write_dot(hpo_plot,file=paste0("hpo_plot_",i))
        
        print(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png saved!"))
        png(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png"), width = 2400, height =1200) 
        plot(hpo_plot)
        dev.off()
        
        
      }else if(nrow(var_allvariant_plot)>=15){
        
        
        
        
        hpo_plot<-onto_plot(hpo, term_sets = list_hpo_patient, label=label_by_term_set,widht = 1,
                            fillcolor = colour_by_frequency,
                            shape = "circle",fontsize = 300,
                            fixedsize = "true")
        
        hpo_plot
        
        #write_dot(hpo_plot,file=paste0("hpo_plot_",i))
        
        print(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png saved!"))
        png(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png"), width = 2400, height =1200) 
        plot(hpo_plot)
        dev.off()
        
        
      }else{
        
        
        hpo_plot<-onto_plot(hpo, term_sets = list_hpo_patient, label=label_by_term_set,widht = 1,
                            fillcolor = colour_by_frequency,
                            shape = "circle",fontsize = 20,
                            fixedsize = "true")
        
        hpo_plot
        
        #write_dot(hpo_plot,file=paste0("hpo_plot_",i))
        
        print(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png saved!"))
        png(paste0("hpo_plot_freq_",sign_gene_test,"_cases.png"), width = 2400, height =1200) 
        plot(hpo_plot)
        dev.off()
        
      }
    }
  }
}else{
      
      setwd(path_to_analysis)
      
      #extract pvalues for disease
      pval_DF_analysis_label<-pval_DF %>% dplyr::filter(analysis.label==disease)
      
      #Show analysis disease label
      sign_analysis_label_read <- as.data.frame(table(pval_DF_analysis_label$analysis))
      sign_analysis_label_read<-as.character(sign_analysis_label_read$Var1)
      print(sign_analysis_label_read)
      
      #create analysis disease folder
      path<-paste0("plots/",sign_analysis_label_read)
      dir.create(file.path(getwd(), path), recursive = TRUE)
      setwd(file.path(getwd(), path))
      getwd()
      
      #Significant pvalues above threshold
      pval_DF_analysis_label_sign<-pval_DF_analysis_label %>% dplyr::filter(p.adjust.fdr<=p.adjust.fdr_threshold)
      
      print(paste0(nrow(pval_DF_analysis_label_sign), " significant gene-disease associations were found in association with ",sign_analysis_label_read))
      
      #read caco definition TSV file
      caco_def<-read_tsv(paste0(path_to_analysis,"/data/",disease,".tsv"))
      
      print(paste0("Reading ",disease,".tsv file"))
      
      #Create loop per each significant gene
      for(sign_gene in unique(pval_DF_analysis_label_sign$gene)){
        
        print(sign_gene)
        
        #Get significant pvalues per that gene
        pval_sign_gene<-pval_DF_analysis_label_sign %>% dplyr::filter(gene==sign_gene)
        
        #Create loop per each significant test
        for(sign_test in pval_sign_gene$test){
          
          sign_gene_test<-paste0(sign_gene,"_",sign_test)
          
          print(sign_gene_test)
          
          #Create folder for each significan gene-test
          dir.create(path = paste0(path_to_analysis,"/plots","/",sign_analysis_label_read,"/",sign_gene_test), recursive = TRUE)
          path<-paste0(path_to_analysis,"/plots","/",sign_analysis_label_read,"/",sign_gene_test)
          setwd(path)
          getwd()
          
          #Extract row for significant gene-test
          pval_sign_gene_test<-pval_sign_gene %>% dplyr::filter(gene==sign_gene) %>% dplyr::filter(test==sign_test)
          
          #Extract all variants in that gene
          gene_all<-exomiser  %>% dplyr::filter(gene==sign_gene)
          
          #Extract caco or caco.denovo info
          if(sign_test=="denovo"){
            
            gene_all$caco<-caco_def$caco.denovo[match(gene_all$sample.id,caco_def$sample.id)]
            
            
          }else{
            
            gene_all$caco<-caco_def$caco[match(gene_all$sample.id,caco_def$sample.id)]
            
          }
          
          
          table(gene_all$caco, useNA = "always")
          
          
          if (file.exists(paste0(path_to_analysis, "/output/probandDataMatrix/probandDataMatrix_", disease, ".tsv"))) {
            
            #####################
            contrib_variants <- read_tsv(paste0(path_to_analysis, "/output/probandDataMatrix/probandDataMatrix_", disease, ".tsv"))
            
            print(paste0("Reading probandDataMatrix_", disease, ".tsv file"))
            
            contrib_variants <- contrib_variants %>%
              dplyr::filter((drop_1 == 1 & drop_2 == 0) | (drop_1 == 0 & drop_2 == 1) | (drop_1 == 0 & drop_2 == 0) | (drop_1 == 0 & is.na(drop_2))) %>%
              dplyr::filter(gene == sign_gene)
            
            # Filter variants based on the test that has been run
            if (sign_test == "LoF") {
              contrib_variants <- contrib_variants %>% dplyr::filter(max.vartype == "LoF")
            } else if (sign_test == "zero80") {
              contrib_variants <- contrib_variants %>% dplyr::filter(max.varscore >= 0.80)
            } else if (sign_test == "denovo") {
              contrib_variants <- contrib_variants %>% dplyr::filter(max.denovo == '1' & fam.structure >= 3)
            } else {
              contrib_variants <- contrib_variants %>% dplyr::filter(max.ccr == 1 & max.varscore >= 0.8)
            }
            
            # Extract caco or caco.denovo info
            if (sign_test == "denovo") {
              contrib_variants$caco <- caco_def$caco.denovo[match(contrib_variants$sample.id, caco_def$sample.id)]
            } else {
              contrib_variants$caco <- caco_def$caco[match(contrib_variants$sample.id, caco_def$sample.id)]
            }
            
            table(contrib_variants$caco, useNA = "always")
            
            cases_contrib <- contrib_variants %>% dplyr::filter(caco == 1)
            controls_contrib <- contrib_variants %>% dplyr::filter(caco == 0)
            
            if (length(unique(cases_contrib$sample.id)) == pval_sign_gene_test$d) {
              print(paste0("Match: ", length(unique(cases_contrib$sample.id)), " variants found == d (", pval_sign_gene_test$d, ")"))
            } else {
              print(paste0("Mismatch: ", length(unique(cases_contrib$sample.id)), " variants found != d (", pval_sign_gene_test$d, ")"))
              break
            }
            
            if (length(unique(controls_contrib$sample.id)) == pval_sign_gene_test$c) {
              print(paste0("Match: ", length(unique(controls_contrib$sample.id)), " variants found == c (", pval_sign_gene_test$c, ")"))
            } else {
              print(paste0("Mismatch: ", length(unique(controls_contrib$sample.id)), " variants found != c (", pval_sign_gene_test$c, ")"))
              break
            }
            
            contrib_df <- contrib_variants %>%
              dplyr::select(sample.id, variant.lift_1, variant.lift_2) %>%
              reshape2::melt(id.vars = "sample.id") %>%
              dplyr::filter(!is.na(value))
            
          }
          
          #Add readable caco definition
          gene_all$caco_definition<-ifelse(gene_all$caco==1,"case",NA)
          gene_all$caco_definition<-ifelse((is.na(gene_all$caco_definition)==TRUE) & gene_all$caco==0,"control",gene_all$caco_definition)
          gene_all$caco_definition<-ifelse((is.na(gene_all$caco_definition)==TRUE),"excluded",gene_all$caco_definition)
          
          table(gene_all$caco_definition, useNA = "always")
          
          #Filter variants based on the test that has been run
          if (sign_test=="LoF"){
            gene_test<- gene_all  %>% dplyr::filter(vartype=="LoF")
          }else if(sign_test=="zero80"){
            gene_test<- gene_all  %>% dplyr::filter(variant.score>=0.80)
          }else if(sign_test=="denovo"){
            gene_test<- gene_all  %>% dplyr::filter(de.novo=='Y' & fam.structure >=3)
          }else{
            gene_test<- gene_all  %>% dplyr::filter(ccr.flag==1 & variant.score>=0.8)
          }
          
          if (exists("contrib_df")) {
            
            gene_test$check <- ifelse(gene_test$sample.id %in% contrib_df$sample.id & gene_test$variant.lift %in% contrib_df$value, "Contrib", "Not contrib")
            
            gene_test$caco_definition <- ifelse(gene_test$caco_definition == "case" & gene_test$check == "Not contrib", "case_drop_varfilter", gene_test$caco_definition)
            gene_test$caco_definition <- ifelse(gene_test$caco_definition == "control" & gene_test$check == "Not contrib", "control_drop_varfilter", gene_test$caco_definition)
            gene_test$check<-NULL
            
          }
          
          table(gene_test$caco_definition, useNA = "always")
          
          # Extract number of cases contributing to the test
          cases <- gene_test %>% dplyr::filter(caco_definition == "case")
          controls <- gene_test %>% dplyr::filter(caco_definition == "control")
          
          if(length(unique(cases$sample.id))==pval_sign_gene_test$d){
            
            print(paste0("Match: ",length(unique(cases$sample.id))," variants found == d (",pval_sign_gene_test$d,")"))
            
          }else{
            
            print(paste0("Mismatch: ",length(unique(cases$sample.id))," variants found != d (",pval_sign_gene_test$d,")"))
            
          }
          
          
          if (length(unique(controls$sample.id)) == pval_sign_gene_test$c) {
            print(paste0("Match: ", length(unique(controls$sample.id)), " variants found == c (", pval_sign_gene_test$c, ")"))
          } else {
            print(paste0("Mismatch: ", length(unique(controls$sample.id)), " variants found != c (", pval_sign_gene_test$c, ")"))
            break
          }
          
          # Extract not contributing variants based on the test that has been run
          if (sign_test == "LoF") {
            not_contributing <- gene_all %>% dplyr::filter(vartype != "LoF")
          } else if (sign_test == "zero80") {
            not_contributing <- gene_all %>% dplyr::filter(variant.score < 0.80)
          } else if (sign_test == "denovo") {
            not_contributing <- gene_all %>% dplyr::filter(de.novo != 'Y' | fam.structure < 3)
          } else {
            not_contributing <- gene_all %>% dplyr::filter(ccr.flag != 1 | variant.score < 0.8)
          }
          
          not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "case", "not_contributing_case", not_contributing$caco_definition)
          not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "control", "not_contributing_control", not_contributing$caco_definition)
          not_contributing$caco_definition <- ifelse(not_contributing$caco_definition == "excluded", "not_contributing_excluded", not_contributing$caco_definition)
          
          #Merge contributing and not contributing variants
          gene_test<-rbind(gene_test,not_contributing)
          
          #Add additional info
          gene_test <-merge(x=gene_test, y=pval_sign_gene_test[,c("gene",colnames(pval_sign_gene_test)[!colnames(pval_sign_gene_test) %in% colnames(gene_test)])], by=c("gene"),all.x=TRUE)
          
          #Order columns
          gene_test<-gene_test %>% dplyr::select(caco_definition,caco,sample.id,hpo.terms,hpo.ids,gene,fam.structure,variant:ccr.flag,vartype,everything()) %>% dplyr::arrange(caco_definition,variant)
          
          #Eventually do some data cleaning
          gene_test$hgvs<-sub("\\|$","",gene_test$hgvs)
          gene_test$exomiser.result.count<-sub("\\|$","",gene_test$exomiser.result.count)
          gene_test$genotypes<-sub("\\|$","",gene_test$genotypes)
          
          
          # Order table by case/control definition
          caco_def_order <- c("case", "case_drop_varfilter", "excluded", "control", "control_drop_varfilter", 
                              "not_contributing_case", "not_contributing_excluded", "not_contributing_control")
          
          gene_test <- gene_test[order(factor(gene_test$caco_definition, levels = unique(caco_def_order))), ]
          
          #Write variants table
          write.table(gene_test,paste0(sign_gene_test,"_",disease,"_variants.tsv"), sep="\t", row.names = FALSE)
          
          print(paste0(sign_gene_test,"_",disease,"_variants.tsv written!"))
          
        }
      }
    }

rm(list = ls())
