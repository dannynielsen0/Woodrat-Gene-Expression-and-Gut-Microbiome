##Maaslin

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(phyloseq)
library(ggplot2)
library(DESeq2)
library(data.table)
BiocManager::install("Maaslin2")
library(Maaslin2)


load("physeq_C.rdata")
physeq_C  #phyloseq object

#remove OTUs with no phylum level resolution
physeq_C <- subset_taxa(physeq_C, Phylum != "NA", check.names=FALSE)



#load DE species by diet interaction caecum genes
DE_interaction_caecum_genes <- read.table("caecum_DE_interaction_result.txt", 
                                          header= TRUE, sep = "\t", row.names = 1)

#load the gene count DF
load("caecum_geneCounts_df.rdata")
count_df_caecum <- t(count_df_caecum)

#subset the counts file to only thse that are in our DE list above
DE_SxD_caecum <- subset(count_df_caecum, rownames(count_df_caecum) %in% 
                          rownames(DE_interaction_caecum_genes))


DE_SxD_caecum <- data.frame(t(DE_SxD_caecum))


#read in the MA max dose, water, and minutes/day data
MA_water <- read.table("MA_maxWater_txt", header=TRUE)
MA_dose <- read.table("MA_maxDose.txt", header = TRUE)
MA_minutes <- read.table("MA_maxMinutes_txt", header=TRUE)

#set up DFs for maaslin
#metadata
meta <- data.frame(physeq_C@sam_data)

#match the phenotypes from above to this metadata file
meta$MA_maxWater <- MA_water$MA_water[match(meta$WR_ID, MA_water$Woodrat_id)]
meta$MA_maxDose <- MA_dose$MA_dose[match(meta$WR_ID, MA_dose$Woodrat_id)]
meta$MA_maxMinutes <- MA_minutes$MA_minutes[match(meta$WR_ID, MA_dose$Woodrat_id)]

#create speciesxdiet variable
meta$species_diet <- paste(meta$Species, meta$Diet_treatment, sep = "_")
meta <- meta[order(meta$WR_ID),]

#save metadata to storage
save(meta, file="Maaslin_meta_phenotypes_Cgenes.rdata")



#create input count data and get the sample ID rownames to match the metadata
input_data <- DE_SxD_caecum
input_data$samID <- rownames(meta) #add rownames to the data
rownames(input_data) <- input_data$samID #this renames the rows to the sample ID that matches the meta


#reorder so it matches the metadata
input_data$samID <- rownames(input_data)
input_data <- input_data[order(input_data$samID),]

#then remove the samID column
input_data <- input_data[,1:367]


#save the input_data for loading in the WGCNA as we will use these
save(input_data, file="SxD_DE_caecum_genes.rdata")

###
#Now have input data and metadata for maaslin
###


#create lepida and bryanti only dfs
meta_lep_FRCA <- subset(meta, meta$Species=="N. lepida" & meta$Diet_treatment=="FRCA")
meta_bry <- subset(meta, meta$Species=="N. bryanti")


#run Maaslin
#for lepida
fit_data = Maaslin2(
  input_data = input_data, 
  input_metadata = meta_lep_FRCA, 
  output = "Maaslin_output/lepida_FRCA", 
  fixed_effects = c("MA_maxWater", "MA_maxDose", "MA_maxMinutes"))

#for bryanti
fit_data = Maaslin2(
  input_data = input_data, 
  input_metadata = meta_bry, 
  output = "Maaslin_output/bryanti", 
  fixed_effects = c("Diet_treatment", "MA_maxWater", "MA_maxDose", "MA_maxMinutes"))






