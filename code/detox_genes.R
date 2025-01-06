#this will summarize the manually annotated detox genes across tissues types 
rm(list=ls())
setwd("/Volumes/MatocqLab/Danny/google_drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#save and load workspace
save.image(file="detox_genes.RData")
load("detox_genes.RData")

#Libraries
# Load the dplyr package
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(png)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(viridis)
  library(wesanderson)
  library(ggtext)


#load in trial metadata
  metadata <- readRDS("lab_trial_final_metadata.RDS")
  
#read in the dataframe of detox genes and numbers of each within lepida and bryanti
  detox_genes <- read.csv("detox_genes.csv", header=TRUE, colClasses = c("character", "character", "integer", "integer"))
  detox_genes$detox_phase <- factor(detox_genes$detox_phase, levels = c("I","II","III"))
  detox_genes <- detox_genes[order(detox_genes$detox_phase),]

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
  all_lep_tissues_deseq <- DESeq(all_lep_tissues)

#pull out the counts for each
  bry_all_counts <- all_bry_tissues@assays@data$counts
  bry_all_counts <- counts(all_bry_tissues_deseq, normalized=TRUE)
  
  lep_all_counts <- all_lep_tissues@assays@data$counts
  lep_all_counts <- counts(all_lep_tissues_deseq, normalized=TRUE)

#using grep to filter count data to only our detox genes of interest

  detox_pattern <- paste(detox_genes$Gene, collapse="|") #create a list of our gene names to grep out below
  detox_pattern <- paste0(detox_pattern, "|TST") #add TST gene to the list
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
  bry_detox_genes_consolidated <- aggregate(. ~ Consolidated_Group, data=bry_detox_genes, sum)
  bry_detox_genes_consolidated <- setNames(data.frame(t(bry_detox_genes_consolidated[,-1])), bry_detox_genes_consolidated[,1]) #transpose

#add a column names aligned and put lep in it
  bry_detox_genes_consolidated <- bry_detox_genes_consolidated %>%
    mutate(aligned = "bry") %>%
    select(aligned, everything())
  
    
#create a conslidated df for lepida-aligned data
#assign names into this new consolidated column
  for (pattern in detox_list) {
    matching_genes <- grepl(pattern, rownames(lep_detox_genes), ignore.case = TRUE)
    lep_detox_genes$Consolidated_Group[matching_genes] <- pattern
  }
  
#conslidate the df on the gene subfamily column  
  lep_detox_genes_consolidated <- aggregate(. ~ Consolidated_Group, data=lep_detox_genes, sum)
  lep_detox_genes_consolidated <- setNames(data.frame(t(lep_detox_genes_consolidated[,-1])), lep_detox_genes_consolidated[,1]) #transpose
  
#add a column names aligned and put lep in it
  lep_detox_genes_consolidated <- lep_detox_genes_consolidated %>%
    mutate(aligned = "lep") %>%
    select(aligned, everything())
  

#Merge the two dataframes 
  combined_detox <- rbind(bry_detox_genes_consolidated,lep_detox_genes_consolidated)
  combined_detox <- bry_detox_genes_consolidated #do for just the bry aligned data
  
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
  
  
#Melt the data frame to long format
  melted_data <- reshape2::melt(combined_detox)
  melted_data <- subset(melted_data, melted_data$aligned=="bry") #subset to only have the bry aligned reads
  
#add phase info for each gene 
  melted_data$phase <- detox_genes$detox_phase[match(melted_data$variable,detox_genes$Gene)]
    
#for loop to create a list that store dataframe for each gene name
  #first, reorder variables for plots and such
  melted_data$phase <- factor(melted_data$phase, levels=c("I","II","III"))
  melted_data$aligned <- factor(melted_data$aligned, levels=c("lep", "bry")) 
  melted_data$Tissue <- factor(melted_data$Tissue, levels=c("FG","SI", "C", "L"))
  melted_data$Species <- factor(melted_data$Species, levels=c("N. lepida", "N. bryanti"))
  melted_data$diet_treatment <- factor(melted_data$diet_treatment, levels=c("PRFA", "FRCA"))
  melted_data$treatment_group <- as.factor(paste(melted_data$Species, melted_data$diet_treatment, sep="_"))
  melted_data$value <- as.numeric(melted_data$value)

 
#add treatmenttype variable that just has specialist or generalist and home and away
  melted_data$treatmenttype <- melted_df$treatmenttype[match(melted_data$treatment_group, melted_df$C_response)]
    

  levels(melted_data$Tissue)
  levels(melted_data$Tissue) <- c("Foregut","Small Intestine", "Caecum", "Liver")
  levels(melted_data$treatment_group)
  levels(melted_data$treatment_group) <- c("Generalist-P","Generalist-NP", "Specialist-NP", "Specialist-P")
  
  melted_data <- melted_data %>%
          group_by(Tissue, variable, phase, treatment_group) %>%
          summarise(mean_value = mean(value))

#generate heatmap with these data
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  modified_palette <- pal 
  modified_palette[1] <- "white"
  
  #change the labels to include 'Phase'
  custom_labels <- function(variable_value) {
    return(dplyr::recode(variable_value,
                  'I' = "Phase I",
                  'II' = "Phase II",
                  'III' = "Phase III"))}
  min_value <- min(melted_data$mean_value, na.rm = TRUE)
  max_value <- max(melted_data$mean_value, na.rm = TRUE)
  
  
#make borders around each facet so the white is shown more, or change the white to grey?
  heatmap <- ggplot(melted_data, aes(y = treatment_group, x = variable, fill = sqrt(mean_value))) +
    geom_tile() + theme_pubclean() + facet_grid(Tissue~phase, scales="free", space="free", labeller = as_labeller(custom_labels)) +
    scale_fill_gradientn(
      colours = modified_palette,
      name = "Square Root\nAverage\nRead Count",
      #na.value = "white",  # Extend color bar to include white for zero
      #guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1)  # Black border around color bar
    ) +
    ggtitle("Detoxification Genes") + theme(plot.title = element_text(hjust = 0.5)) +
    #scale_fill_gradientn(colors=heat_colors, na.value = NA, name = "Square Root(Average Read Count)") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1),
          legend.position = "right", legend.direction = "vertical") +
    # scale_y_discrete(labels = function(x) {
    #   sapply(x, function(label) {
    #     parts <- strsplit(as.character(label), "-", fixed = TRUE)[[1]]
    #     if (length(parts) > 1) {
    #       parse(text = sprintf("italic(%s) ~ '-' ~ %s", sQuote(parts[1], FALSE), sQuote(parts[2], FALSE)))
    #     } else {
    #       label
    #     }
    #   })
    # }) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 2)) +
    labs(x = "", y = "", fill = "Value") +
    theme(text = element_text(size = 30),
          axis.title = element_text(size = 30),
          axis.text.y = element_text(),  # Ensure y-axis text is formatted to show expressions
          axis.text = element_text(size = 30),
          strip.text = element_text(size = 30, face = "bold"))  # Increase facet heading sizes and make them bold
  
  heatmap <- ggplot(melted_data, aes(x = reorder(variable, as.numeric(factor(phase))), y = factor(Tissue, levels=rev(levels(Tissue))), fill = sqrt(mean_value))) +
     geom_tile() + theme_pubclean() + facet_grid(sp_diet~phase, scales="free", space="free", labeller = as_labeller(custom_labels)) +
    scale_fill_gradientn(colours = modified_palette,name = "Square Root\n(Average Read Count)") + 
    #scale_fill_gradientn(colors=heat_colors, na.value = NA, name = "Square Root(Average Read Count)") +
    theme(axis.text.x = element_text(angle = 75, hjust = 1), 
          legend.position = "right", legend.direction = "vertical") +
    labs(x = "", y = "", fill = "Value") +
     theme(text = element_text(size = 30),
           axis.title = element_text(size = 30),
           axis.text = element_text(size = 30),
           strip.text = element_text(size = 30, face = "bold")) # Increase facet heading sizes and make them bold
  
  ggsave(plot=heatmap, "../Lab-diet-trial-16S-analysis/figures/detox_genes_sqrt_count_heatmap_v2.pdf", width = 30, height =20 , device='pdf', dpi=700)
  
  
  # Plot the data with a loop through each gene family
  options(scipen=1000)  
  
  ggplot_list <- list()
  
  # Iterate through gene name
  for (gene in detox_list) {
    # Subset data for each gene
    gene_subset <- melted_data[melted_data$variable == gene, ]
    
    # Check if the subset is not empty
    if (nrow(gene_subset) > 0) {
      # Calculate group averages
      #group_averages <- aggregate(value ~ Tissue + Species + diet_treatment, data = gene_subset, FUN = mean)
      
      # Create ggplot for the current gene
      gg <- ggplot(gene_subset, aes(x = interaction(WR_ID,Species, diet_treatment), y = value, fill=aligned)) +
        # Add dashed line for group average in each facet
        #geom_hline(data = group_averages, aes(yintercept = value, linetype = "Group Average"), 
         #          color = "black", linetype = "dashed", linewidth = 0.4) +
        geom_bar(stat = "identity") + #guides(fill=guide_legend(title="Aligned")) +
        facet_grid(Species ~ Tissue + diet_treatment, scales="free_x", space = "free_x") +
        labs(title = paste("Stacked Bar Plot for", gene, "Phase-", unique(gene_subset$phase)), x = "", y = "Total Reads") +
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

  
  
  
