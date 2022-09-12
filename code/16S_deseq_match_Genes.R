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
library(limma)
library(factoextra)

#also need to install btools for estimating diversity
devtools::install_github('twbattaglia/btools')
library(btools)


#load in phyloseq object for Caecum data

physeq_C <- readRDS("physeq_C.rds")


#reorder factors
physeq_C@sam_data$Species <- factor(physeq_C@sam_data$Species, levels=c("N. lepida", "N. bryanti"))

#make sp_diet group
physeq_C@sam_data$sp_diet <- paste(physeq_C@sam_data$Species, physeq_C@sam_data$Diet_treatment, sep="_")

#glom at genus level

physeq_C_family <- tax_glom(physeq_C, taxrank = "Family")

#change rownames to microbial family name
tax_Df <- data.frame(physeq_C_family@tax_table)
rownames(physeq_C_family@otu_table) <- tax_Df$Family[match(rownames(physeq_C_family@otu_table), rownames(tax_Df))] #get the genus level microbial taxonomy


#construct the deseq object
dds = phyloseq_to_deseq2(physeq_C_family, ~0 + sp_diet)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)



# perform a PCA on the data in assay(x) for the selected genes
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)




biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 25), #draw top 25 arrows
                          addEllipses = TRUE,
                          habillage = dds$sp_diet,
                          col.ind = dds$Species,
                          ellipse.level=0.95,
                          palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                          geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                          ind.shape = dds$Diet,
                          ind.fill = dds$Species,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Microbiome of Woodrat Caecum") +
      theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

ggsave(plot=biplot, "../Lab-diet-trial-16S-analysis/figures/microbiome_PCA_familyLevel.jpg",
       width = 10, height =10 , device='jpg', dpi=500)


