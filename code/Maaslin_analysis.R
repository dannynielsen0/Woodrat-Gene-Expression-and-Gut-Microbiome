##Maaslin analysis for lab diet experiments to explore phenotypes to genes in caecum

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(phyloseq)
library(ggplot2)
library(DESeq2)
library(data.table)
BiocManager::install("Maaslin2")
library(Maaslin2)

#load input data 
load("SxD_DE_caecum_genes.rdata")

#load metadata
load("Maaslin_meta_phenotypes_Cgenes.RData")

#run Maaslin on lepida and bryanti separately

#create lepida and bryanti only dfs
meta_lep <- subset(meta, meta$Species=="N. lepida")

meta_bry <- subset(meta, meta$Species=="N. bryanti")



#for lepida
fit_data = Maaslin2(
  input_data = input_data, 
  input_metadata = meta_lep, 
  output = "Maaslin_output/lepida", 
  fixed_effects = c("Diet_treatment", "MA_maxWater", "MA_maxDose", "MA_maxMinutes"))

#for bryanti
fit_data = Maaslin2(
  input_data = input_data, 
  input_metadata = meta_bry, 
  output = "Maaslin_output/bryanti", 
  fixed_effects = c("Diet_treatment", "MA_maxWater", "MA_maxDose", "MA_maxMinutes"))






