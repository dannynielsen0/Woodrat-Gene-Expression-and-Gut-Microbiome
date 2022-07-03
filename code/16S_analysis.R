
###Microbiome analysis for chapter 3 - feeding trials

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for such objects
library("biomformat") # perhaps unnecessary 
library(ggplot2)
library(picante) 
library(decontam)
library(WGCNA)
library(DESeq2)
library(data.table)
library("EnvStats") #plotting sample size
library("ggpubr") #plotting stats
library(vegan)


#also need to install btools for estimating diversity
devtools::install_github('twbattaglia/btools')
library(btools)


#load in phyloseq object for Caecum data

load("physeq_C.rdata")
physeq_C  #phyloseq object

#remove potential contaminants
physeq_C <- subset_taxa(physeq_C, Family != "Mitochondria" | Class != "Chloroplast")



#check sequencing depth
sdt = data.table(as(sample_data(physeq_C), "data.frame"),
                 TotalReads = sample_sums(physeq_C), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth



#calculate Faith's phylogenetic diversity

physeq_C@sam_data$PD <- estimate_pd(physeq_C) #uses the btools package

#Then plot 
PD_div_box_plot <- ggplot(data=physeq_C@sam_data, aes(x=physeq_C@sam_data$Diet_treatment,y=physeq_C@sam_data$PD$PD)) +
  geom_boxplot() + facet_wrap(~physeq_C@sam_data$Species) + geom_jitter() + stat_n_text() + theme_bw() +
  geom_signif(test = "wilcox.test", y_position = 13, map_signif_level = TRUE, comparisons = list(c("PRFA", "FRCA"))) + ylab("Diversity") + xlab("Diet Treatment")

PD_div_box_plot


#convert to RRA 

#This will convert reads counts to RRA
physeq_RRA <- transform_sample_counts(physeq_C, function(x) x/sum(x)) #converts to relative abundance

#check it now
physeq_RRA@otu_table[1:3]

#construct ordinations for analysis and plotting

PCoA_bray <- ordinate(physeq_RRA, method = "PCoA", distance = "bray") #ordination using bray-curtis distances
PCoA_wunifrac <- ordinate(physeq_RRA, method = "PCoA", distance = "wunifrac") #ordination using weighted unifrac
PCoA_unifrac <- ordinate(physeq_RRA, method = "PCoA", distance = "unifrac") #unweighted unifrac

#plot bray
PCoA_plot_bray <- plot_ordination(physeq_RRA, PCoA_bray, color = "Diet_treatment", axes = 1:2)

plot_bray <- PCoA_plot_bray + geom_point(size = 5) +
  scale_color_manual(values= c("forestgreen", "maroon")) + facet_wrap(~Species) +
  scale_shape_manual(values=c(17,15)) + theme_bw() 

print(plot_bray)

#plot weighted unifrac

PCoA_plot_wunifrac <- plot_ordination(physeq_RRA, PCoA_wunifrac, color = "Diet_treatment", axes = 1:2)

plot_wunifrac <- PCoA_plot_wunifrac + geom_point(size = 5) +
  scale_color_manual(values= c("forestgreen", "maroon")) + facet_wrap(~Species) +
  scale_shape_manual(values=c(17,15)) + theme_bw() 

print(plot_wunifrac)


#plot unweighted unifrac

PCoA_plot_unifrac <- plot_ordination(physeq_RRA, PCoA_unifrac, color = "Diet_treatment", axes = 1:2)

plot_unifrac <- PCoA_plot_unifrac + geom_point(size = 5) +
  scale_color_manual(values= c("forestgreen", "maroon")) + facet_wrap(~Species) +
  scale_shape_manual(values=c(17,15)) + theme_bw() 

print(plot_unifrac)


#perform PERMANOVA

perm_metadata <- as(sample_data(physeq_RRA), "data.frame")

#then we calculate a distance matrix - do bray-curtis, and both weighted and unweighted unifrac
distmat_bray <- phyloseq::distance(physeq_RRA, method="bray") #use bray here. note, the :: tells R to use the phyloseq distance function
distmat_unifrac <- phyloseq::distance(physeq_RRA, method="unifrac") #use bray here. note, the :: tells R to use the phyloseq distance function
distmat_wunifrac <- phyloseq::distance(physeq_RRA, method="wunifrac") #use bray here. note, the :: tells R to use the phyloseq distance function

#run the analyses
permanova_bray <- adonis2(distmat_bray ~ Species + Diet_treatment + Sex + Species*Diet_treatment, data=perm_metadata, permutations = 999)
permanova_unifrac <- adonis2(distmat_unifrac ~ Species + Diet_treatment + Sex + Species*Diet_treatment, data=perm_metadata, permutations = 999)
permanova_wunifrac <- adonis2(distmat_wunifrac ~ Species + Diet_treatment + Sex + Species*Diet_treatment, data=perm_metadata, permutations = 999)

permanova_bray
permanova_unifrac
permanova_wunifrac

