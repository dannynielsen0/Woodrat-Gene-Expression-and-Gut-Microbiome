#This will load Qiime data and create phyloseq object in R environment for 16S analysis

#load libraries
library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for 16S data
library(ggplot2) #plotting
library(picante) #diversity estimates
library(decontam) #remove contaminants based on those present in Blank

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

# import QIIME2 data to phyloseq

#need qiime2R package to convert Qiime files to phyloseq
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

physeq <- qza_to_phyloseq(features="table.qza",tree="fasttree-tree-rooted.qza",
                          taxonomy = "taxonomy.qza", metadata = "metadata.tsv")

