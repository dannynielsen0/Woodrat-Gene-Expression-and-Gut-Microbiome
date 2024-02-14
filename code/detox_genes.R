
#Libraries
# Load the dplyr package
  library(dplyr)
  library(tibble)


#this will summarize the manually annotated detox genes across tissues types 
  rm(list=ls())
  setwd("~/Library/CloudStorage/GoogleDrive-dannynielsen@nevada.unr.edu/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#read in the dataframe of detox genes and numbers of each within lepida and bryanti
  detox_genes <- read.csv("detox_genes.csv", header=TRUE, colClasses = c("character", "integer", "integer"))

#read in the count data for bry and secondary lep. aligned data
#these data are in bryFirst_counts_files and lepSecond_counts_files

#From the bryanti first pass alighments and counts, choose which tissue to analyze, change bry to lep for lepida aligned data
  all_bry_directory <- "bryFirst_counts_files/all_bry_counts"


#From the lepida second pass alighments and counts, choose which tissue to analyze, change bry to lep for lepida aligned data
  all_lep_directory <- "lepSecond_counts_files/all_lep_counts"

#get  names of htseq-count output files
  bry_files <-  list.files(all_bry_directory, pattern = ".tsv")
  lep_files <- list.files(all_lep_directory, pattern = ".tsv")

#load metadata
  sampleTable <- read.csv("Meta_data/woodrat_file_metadata.csv", header = TRUE)

#build group variable to facilitate contrasts (rather than specifying interactions model)
#approach detailed in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
  sampleTable$group <- factor(paste0(sampleTable$Species,"_",sampleTable$Diet))

#build dds object for storing all data, design is a linear model formula given variables to test
#one each for bry and lep aligned counts
  all_bry_tissues <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                                   directory = all_bry_directory,
                                                   design= ~ 0 + group)

  all_lep_tissues <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                              directory = all_lep_directory,
                                              design= ~ 0 + group)


#pull out the counts for each
  bry_all_counts <- all_bry_tissues@assays@data$counts
  lep_all_counts <- all_lep_tissues@assays@data$counts

#using grep to filter count data to only our detox genes of interest

  detox_pattern <- paste(detox_genes$Gene, collapse="|") #create a list of our gene names to grep out below
  remove_pattern <- paste("Fmod", "-like", "pseudo", sep="|") #pattern matches that need to be removed
  
#create a detox only dataset for bryanti
  bry_detox_genes <- as.data.frame(bry_all_counts[grep(detox_pattern, row.names(bry_all_counts), ignore.case=TRUE),]) #use grep to create a df of counts of only detox genes
  bry_detox_genes <- (bry_detox_genes[grep(remove_pattern, row.names(bry_detox_genes), invert = TRUE, ignore.case=TRUE),]) #use grep to create a df of counts of only detox genes
  
#create a detox only dataset for lepida
  lep_detox_genes <- as.data.frame(lep_all_counts[grep(detox_pattern, row.names(lep_all_counts), ignore.case = TRUE),])
  lep_detox_genes <- (lep_detox_genes[grep(remove_pattern, row.names(lep_detox_genes), invert = TRUE, ignore.case = TRUE),])

#clean up gene names in rows (i.e. Nbry, Nlep, etc...)
  rownames(bry_detox_genes) <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "", rownames(bry_detox_genes))
  rownames(lep_detox_genes) <- gsub("NLEP_|gene-|_Nmacr|_Nlep", "", rownames(lep_detox_genes))
  

#unlist the detox_patterns
  detox_list <- unlist(strsplit(detox_pattern, "|", fixed = TRUE))

#create a conslidated df for bryanti-aligned data
bry_detox_genes$Consolidated_Group <- NA #create an empty column for the consolidated names
#assign gene subfamily names into new consolidated column  
  for (pattern in detox_list) {
    matching_genes <- grepl(pattern, rownames(bry_detox_genes), ignore.case = TRUE)
    bry_detox_genes$Consolidated_Group[matching_genes] <- pattern
  }

#conslidate the df on the gene subfamily column  
  bry_detox_genes_conslidated <- aggregate(. ~ Consolidated_Group, data=bry_detox_genes, sum)
  
  
#create a conslidated df for lepida-aligned data
#assign names into this new consolidated column
  for (pattern in detox_list) {
    matching_genes <- grepl(pattern, rownames(lep_detox_genes), ignore.case = TRUE)
    lep_detox_genes$Consolidated_Group[matching_genes] <- pattern
  }
  
#conslidate the df on the gene subfamily column  
  lep_detox_genes_conslidated <- aggregate(. ~ Consolidated_Group, data=lep_detox_genes, sum)
  

  
     
###Everything below is notes from other code to pull from for above
  #to look for gene families
  View(all_counts[grep("NBRY_CYP[1-3]([^0-9]|$)", row.names(all_counts)),]) #replace SULT with desired gene
View(caecum.int.DEGs[grep("SULT", row.names(caecum.int.DEGs)),])
View(liver_spec_DEGs[grep("NBRY_GST", row.names(liver_spec_DEGs)),])

all_cyps <- t(all_counts[grep("_CYP|-cyp", row.names(all_counts)),])
all_sults <- t(all_counts[grep("_SULT", row.names(all_counts)),])
all_Ugts <- t(all_counts[grep("-Ugt", row.names(all_counts)),])
all_Ugts <- cbind(all_Ugts,sampleTable)
all_cyps <- cbind(all_cyps,sampleTable) #add sample info
all_cyps <- all_cyps[,c(121:129,1:120)]
all_Ugts <- all_Ugts[,c(30:38,1:29)]
write.csv(all_cyps_w, "cyps.csv")
#save all_counts object
saveRDS(all_counts, file = "all_gene_counts.RDS") 




