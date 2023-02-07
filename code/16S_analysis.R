
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
library(qiime2R)


#also need to install btools for estimating diversity
devtools::install_github('twbattaglia/btools')
library(btools)


#load in phyloseq object for Caecum data

physeq <- qza_to_phyloseq(features="qiime_output/table.qza",tree="qiime_output/fasttree-tree-rooted.qza",
                          taxonomy = "qiime_output/taxonomy.qza", metadata = "qiime_output/C_FG_metadata.txt")

#change variables to factors
physeq@sam_data[] <- lapply(physeq@sam_data, factor)


sum(sample_sums(physeq@otu_table)) #635,508



#remove potential contaminants
physeq<- subset_taxa(physeq, Family != "Mitochondria" & Order != "Chloroplast" & Kingdom != "d__Eukaryota")

#remove singletons
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq)

#make group species_diet 
physeq@sam_data$sp_diet <- paste(physeq@sam_data$Species, physeq@sam_data$diet_treatment, sep="_")

#rename RHCA to FRCA
levels(physeq@sam_data$diet_treatment) <- c("PRFA", "FRCA")

#make separate for C and FG samples
physeq_C <- subset_samples(physeq, Gut_region=="Caecum")
physeq_FG <- subset_samples(physeq, Gut_region=="Foregut")

#check sequencing depth
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth


#calculate Faith's phylogenetic diversity

physeq_C@sam_data$PD <- estimate_pd(physeq_C) #uses the btools package
physeq_FG@sam_data$PD <- estimate_pd(physeq_FG)#uses the btools package


#Then plot using wilcox test for differences in diversity between groups - plot both PD & richness
PD_div_box_plot <- ggplot(data=physeq_C@sam_data, aes(x=physeq_C@sam_data$diet_treatment,y=physeq_C@sam_data$PD$SR)) +
  geom_boxplot() + facet_wrap(~physeq_C@sam_data$Species) + geom_jitter() + stat_n_text() + theme_bw() +
  geom_signif(test = "wilcox.test", y_position = 275, map_signif_level = TRUE, comparisons = list(c("PRFA", "FRCA"))) +
  ylab("Microbial Richness of Caecum") + xlab("Diet Treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 24)) + 
  theme(legend.text = element_text(size = 20, face = "italic")) +
  theme(legend.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20)) + theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20))


PD_div_box_plot

#save plot
ggsave(plot=PD_div_box_plot, "../Lab-diet-trial-16S-analysis/figures/microbial_richness_C.jpg", width = 10, height =8 , device='jpg', dpi=500)


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
PCoA_plot_bray <- plot_ordination(physeq_RRA, PCoA_bray, color = "Species", shape = "Diet_treatment", axes = 1:2)

plot_bray <- PCoA_plot_bray + geom_point(size = 5) +
  scale_color_manual(values= c("maroon","forestgreen")) + #facet_wrap(~Species) +
  scale_shape_manual(values=c(16,2)) + theme_bw()+ ggtitle("Bray-Curtis") +
  theme(plot.title = element_text(hjust = 0.5, size = 24)) + 
  theme(legend.text = element_text(size = 20, face = "italic")) +
  theme(legend.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20)) + theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20))

print(plot_bray)

#save plot
ggsave(plot=plot_bray, "../Lab-diet-trial-16S-analysis/figures/bray_PCoA.jpg", width = 10, height =8 , device='jpg', dpi=500)


#plot weighted unifrac

PCoA_plot_wunifrac <- plot_ordination(physeq_RRA, PCoA_wunifrac, color = "Species", shape = "Diet_treatment", axes = c(1,2))

plot_wunifrac <- PCoA_plot_wunifrac + geom_point(size = 5) +
  scale_color_manual(values= c("maroon","forestgreen")) + #facet_wrap(~Species) +
  scale_shape_manual(values=c(16,2)) + theme_bw() + ggtitle("Weighted UniFrac") +
  theme(plot.title = element_text(hjust = 0.5, size = 24)) + 
  theme(legend.text = element_text(size = 20, face = "italic")) +
  theme(legend.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20)) + theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20))

print(plot_wunifrac)
#save plot
ggsave(plot=plot_wunifrac, "../Lab-diet-trial-16S-analysis/figures/weighted_unifrac_PCoA.jpg", width = 10, height =8 , device='jpg', dpi=500)


#plot unweighted unifrac

PCoA_plot_unifrac <- plot_ordination(physeq_RRA, PCoA_unifrac, color = "Species", shape = "Diet_treatment", axes = 1:2)

plot_unifrac <- PCoA_plot_unifrac + geom_point(size = 5) +
  scale_color_manual(values= c("maroon","forestgreen")) + #facet_wrap(~Species) +
  scale_shape_manual(values=c(16,2)) + theme_bw() + ggtitle("Unweighted UniFrac") +
  theme(plot.title = element_text(hjust = 0.5, size = 24)) + 
  theme(legend.text = element_text(size = 20, face = "italic")) +
  theme(legend.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20)) + theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20))
  

print(plot_unifrac)
#save plot
ggsave(plot=plot_unifrac, "../Lab-diet-trial-16S-analysis/figures/unifrac_PCoA.jpg", width = 10, height =8 , device='jpg', dpi=500)


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
write.table(permanova_bray, "../Lab-diet-trial-16S-analysis/figures/bray_permanova.txt", sep= "\t")

permanova_unifrac
write.table(permanova_unifrac, "../Lab-diet-trial-16S-analysis/figures/unifrac_permanova.txt", sep= "\t")

permanova_wunifrac
write.table(permanova_wunifrac, "../Lab-diet-trial-16S-analysis/figures/weighted_unifrac_permanova.txt", sep= "\t")

