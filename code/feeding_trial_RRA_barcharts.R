##RRA barchart plotting for woodrat lab geeding experiment caecum 16S data

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(RColorBrewer)
library(phyloseq)     # filtering and utilities for such objects

physeq <- qza_to_phyloseq(features="qiime_output/table.qza",tree="qiime_output/fasttree-tree-rooted.qza",
                          taxonomy = "qiime_output/taxonomy.qza", metadata = "qiime_output/C_FG_metadata.txt")

sum(sample_sums(physeq_FG@otu_table)) #635,508


#remove potential contaminants
physeq<- subset_taxa(physeq, Family != "Mitochondria" & Order != "Chloroplast" & Kingdom != "d__Eukaryota")

#remove singletons
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq)

#make group species_diet 
physeq@sam_data$sp_diet <- paste(physeq@sam_data$Species, physeq@sam_data$diet_treatment, sep="_")



#make separate for C and FG samples
physeq_C <- subset_samples(physeq, Gut_region=="Caecum")
physeq_FG <- subset_samples(physeq, Gut_region=="Foregut")


#organize data and plot RRA barcharts

########To plot the desired taxanomic level, change each level between this line and the next string of #s below

y1 <- tax_glom(physeq_FG, taxrank = 'Genus') # agglomerate taxa at desired level
y2 = suppressWarnings(merge_samples(y1, "sp_diet")) # merge samples on sample variable of interest; suppress warnings is to prevent from failing
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- suppressWarnings(psmelt(y3)) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Genus[y4$Abundance < 0.05] <- "Genera < 5% abund." #rename genera with < .25% abundance

#split back out the location and staph into unique columns
y4 <- cbind(y4, read.table(text=y4$Sample, sep="_", header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
y4$Genus <- factor(y4$Genus, levels=rev(unique(y4$Genus)))

#change order so lepida on left
y4$col1 <- factor(y4$col1, levels = c("N. lepida", "N. bryanti"))


#set color palette to accommodate the number of genera
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(y4$Genus)))


#plot
p <- ggplot(data=y4, aes(x=col2, y=Abundance, fill=Genus, Genus =Genus)) +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill=guide_legend(nrow=5)) + 
  scale_x_discrete(limits=rev(levels(as.factor(y4$col2)))) + facet_wrap(~col1) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  guides(fill = guide_legend(reverse = TRUE)) + ggtitle("Foregut") +
  ylab("16S relative read abundance") + xlab("")
p


#Genus
ggsave(plot=p, "../Lab-diet-trial-16S-analysis/figures/foregut_phyla_RRA.jpg", width = 14, height =8 , device='jpg', dpi=500)
########

#This is just a quick plotting of data to dummy check my above resulting figures,
#change names and taxonomic levels as needed

clostridium <- subset_taxa(physeq_C, Phylum=="Bacteroidota")

plot_bar(clostridium) + facet_wrap(~sp_diet)
