#This will load Qiime data and create phyloseq object in R environment for 16S analysis

rm(list=ls())

setwd("~/Desktop/lab_trials_16S/caecum_fastq/qiime_output")

#load libraries
library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for 16S data
library(ggplot2) #plotting
library(picante) #diversity estimates
library(decontam) #remove contaminants based on those present in Blank


# import QIIME2 data to phyloseq

#need qiime2R package to convert Qiime files to phyloseq
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R", force=TRUE)
library(qiime2R)

physeq <- qza_to_phyloseq(features="table.qza",tree="fasttree-tree-rooted.qza",
                          taxonomy = "taxonomy.qza", metadata = "../caecum_metadata.tsv")

#remove BLANK and rename UserRef to Diet_treatment
physeq <- subset_samples(physeq, Sample.Type != "BLANK")
colnames(physeq@sam_data)[8] <- "Diet_treatment"

#remove potential contaminants and singletons

#remove potential contaminants
physeq <- subset_taxa(physeq, Family != "Mitochondria" | Class != "Chloroplast")

#remove singletons
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq)

#save this phyloseq object as physeq_C

saveRDS(physeq, "/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data/physeq_C.rds")


