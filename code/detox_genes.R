
#Libraries
# Load the dplyr package
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)

#this will summarize the manually annotated detox genes across tissues types 
  rm(list=ls())
  setwd("~/Library/CloudStorage/GoogleDrive-dannynielsen@nevada.unr.edu/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#load in trial metadata
  metadata <- readRDS("lab_trial_final_metadata.RDS")
  
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

#run deseq
  all_bry_tissues_deseq <- DESeq(all_bry_tissues)

#pull out the counts for each
  bry_all_counts <- all_bry_tissues@assays@data$counts
  lep_all_counts <- all_lep_tissues@assays@data$counts

#using grep to filter count data to only our detox genes of interest

  detox_pattern <- paste(detox_genes$Gene, collapse="|") #create a list of our gene names to grep out below
  remove_pattern <- paste("Fmod", "-like", "pseudo", sep="|") #pattern matches that need to be removed
  
#create a detox only dataset for bryanti
  bry_detox_genes <- as.data.frame(bry_all_counts[grep(detox_pattern, row.names(bry_all_counts), ignore.case=TRUE),]) #use grep to create a df of counts of only detox genes
  bry_detox_genes <- (bry_detox_genes[grep(remove_pattern, rownames(bry_detox_genes), invert = TRUE, ignore.case=TRUE),]) #use grep to create a df of counts of only detox genes
  
#create a detox only dataset for lepida
  lep_detox_genes <- as.data.frame(lep_all_counts[grep(detox_pattern, row.names(lep_all_counts), ignore.case = TRUE),])
  lep_detox_genes <- (lep_detox_genes[grep(remove_pattern, rownames(lep_detox_genes), invert = TRUE, ignore.case = TRUE),])

#clean up gene names in rows (i.e. Nbry, Nlep, etc...)
  rownames(bry_detox_genes) <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "", rownames(bry_detox_genes))
  rownames(lep_detox_genes) <- gsub("NLEP_|gene-|_Nmacr|_Nlep", "", rownames(lep_detox_genes))
  

#unlist the detox_patterns
  detox_list <- unlist(strsplit(detox_pattern, "|", fixed = TRUE))

#create a conslidated df for bryanti-aligned data
#assign gene subfamily names into new consolidated column  
  for (pattern in detox_list) {
    matching_genes <- grepl(pattern, rownames(bry_detox_genes), ignore.case = TRUE)
    bry_detox_genes$Consolidated_Group[matching_genes] <- pattern
  }

#conslidate the df on the gene subfamily column  
  bry_detox_genes_conslidated <- aggregate(. ~ Consolidated_Group, data=bry_detox_genes, sum)
  bry_detox_genes_conslidated <- setNames(data.frame(t(bry_detox_genes_conslidated[,-1])), bry_detox_genes_conslidated[,1]) #transpose

#add a column names aligned and put lep in it
  bry_detox_genes_conslidated <- bry_detox_genes_conslidated %>%
    mutate(aligned = "bry") %>%
    select(aligned, everything())
  
    
#create a conslidated df for lepida-aligned data
#assign names into this new consolidated column
  for (pattern in detox_list) {
    matching_genes <- grepl(pattern, rownames(lep_detox_genes), ignore.case = TRUE)
    lep_detox_genes$Consolidated_Group[matching_genes] <- pattern
  }
  
#conslidate the df on the gene subfamily column  
  lep_detox_genes_conslidated <- aggregate(. ~ Consolidated_Group, data=lep_detox_genes, sum)
  lep_detox_genes_conslidated <- setNames(data.frame(t(lep_detox_genes_conslidated[,-1])), lep_detox_genes_conslidated[,1]) #transpose
  
#add a column names aligned and put lep in it
  lep_detox_genes_conslidated <- lep_detox_genes_conslidated %>%
    mutate(aligned = "lep") %>%
    select(aligned, everything())
  

#Merge the two dataframes 
  combined_detox <- rbind(bry_detox_genes_conslidated,lep_detox_genes_conslidated)
  
#add columns for WR ID and tissue
  combined_detox <- combined_detox %>%
    rownames_to_column(var = "Row_Names") %>%
    separate(Row_Names, into = c("WR_ID", "Tissue"), sep="_", extra="merge", remove=FALSE) %>%
    select(WR_ID, Tissue, everything()) %>%
    select(-Row_Names)

#clean up the tissue column
  combined_detox$Tissue <- sub("_.*$", "", combined_detox$Tissue)
  combined_detox$Tissue <- sub("S\\d{1,3}", "L", combined_detox$Tissue)#change to L for liver
  
#add species, diet treatement, and dose values to the df
  combined_detox$Species <- metadata$Species[match(combined_detox$WR_ID,metadata$WR_ID)]
  combined_detox <- combined_detox[,c(ncol(combined_detox), 1:ncol(combined_detox)-1)]#move the species column to be second
  combined_detox$diet_treatment <- metadata$Diet_treatment[match(combined_detox$WR_ID,metadata$WR_ID)]
  
  
# Melt the data frame to long format
  melted_data <- reshape2::melt(combined_detox)
  
#filter out a specific gene subfamily = ABCA
#for loop to create a list that store dataframe for each gene name
  #first, reorder variables for plots and such
  melted_data$aligned <- factor(melted_data$aligned, levels=c("lep", "bry")) 
  melted_data$Tissue <- factor(melted_data$Tissue, levels=c("FG","SI", "C", "L"))
  melted_data$Species <- factor(melted_data$Species, levels=c("N. lepida", "N. bryanti"))
  melted_data$diet_treatment <- factor(melted_data$diet_treatment, levels=c("PRFA", "FRCA"))
  melted_data$value <- as.numeric(melted_data$value)
  

# Plot the data with a loop through each gene family
  options(scipen=1000)  
  
  ggplot_list <- list()
  
  # Iterate through gene names
  for (gene in detox_list) {
    # Subset data for each gene
    gene_subset <- melted_data[melted_data$variable == gene, ]
    
    # Check if the subset is not empty
    if (nrow(gene_subset) > 0) {
      # Calculate group averages
      group_averages <- aggregate(value ~ Species + Tissue + diet_treatment, data = gene_subset, FUN = mean)
      
      # Create ggplot for the current gene
      gg <- ggplot(gene_subset, aes(x = interaction(WR_ID,Species, diet_treatment), y = value, fill=aligned)) +
        # Add dashed line for group average in each facet
        geom_hline(data = group_averages, aes(yintercept = value, linetype = "Group Average"), 
                   color = "black", linetype = "dashed", linewidth = 0.4) +
        geom_bar(stat = "identity") + #guides(fill=guide_legend(title="Aligned")) +
        facet_grid(Species ~ Tissue + diet_treatment, scales="free_x", space = "free_x") +
        labs(title = paste("Stacked Bar Plot for", gene), x = "", y = "Total Reads") +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank()) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
        
      # Store the ggplot object in the list
      ggplot_list[[gene]] <- gg
    } else {
      # Print a message indicating that the subset is empty
      cat("Skipping gene", gene, "as there are no observations.\n")
    }
  }
  
  # Save all plots to a single PDF file
  pdf("stacked_bar_plots.pdf")
  
  # Loop through the ggplot list and print each plot to the PDF
  for (gene_plot in ggplot_list) {
    print(gene_plot)
  }
  
  # Close the PDF file
  dev.off()

  
  
  
