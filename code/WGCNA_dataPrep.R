
#WGCNA analysis for caecum gene and microbiome 

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(phyloseq)
library(ggplot2)
library(DESeq2)
library(data.table)
library(WGCNA)

#load in phyloseq object for Caecum data

load("physeq_C.rdata")
physeq_C  #phyloseq object

physeq_C_removeNA <- subset_taxa(physeq_C, Phylum != "NA", check.names=FALSE)

#load in the RNAseq gene count data
load("caecum_geneCounts_df.rdata")

#load in the DE genes relateed to interaction with other organisms 
#and the DE microbes from lepida that differed between diets - we'll subset the genes and microbes with these

DE_microbes <- read.table("DESeq_lep_diet.tsv", sep="\t")
DE_genes <- read.table("DE_gene_functions.txt", sep = "\t", header = TRUE)

#subset the phyloseq object to only OTUs in DE_microbes
my_sub <- subset(physeq_C_removeNA@otu_table, rownames(physeq_C_removeNA@otu_table) %in% DE_microbes$x)
DE_microbes_physeq <- merge_phyloseq(my_sub, physeq_C_removeNA@tax_table, physeq_C_removeNA@sam_data)

#subset gene counts for caecum
count_caecum <- data.frame(t(count_df_caecum))
counts_caecum <- subset(tcount_caecum[,1:29], rownames(count_caecum) %in% DE_genes$GENE_ID)
counts_caecum <- t(counts_caecum)

#check the order of samples and order such that is in ascending to match the gene counts df
df <- data.frame(DE_microbes_physeq@sam_data) #convert to df
df <-  df[order(df$WR_ID),] # order by WR_ID
df$sampleID <- rownames(df) # add rownames as variable for matching

#create count df of 16S reads
microbe_dat <- data.frame(t(DE_microbes_physeq@otu_table)) #convert otu table to df and transform
microbe_dat$sampleID <- rownames(microbe_dat) #add sampleID as variable so we can match order below

#match the order here so it matches the gene counts 
microbe_dat <- microbe_dat[order(match(microbe_dat$sampleID,df$sampleID)),]


###We should have the two df's in the same sample order now and we can merge into on giant df

gene_microbe_counts <- cbind(microbe_dat,counts_caecum)
gene_microbe_counts$sampleID <- NULL #drop the sampldeID column now so that we only have count data for columns
