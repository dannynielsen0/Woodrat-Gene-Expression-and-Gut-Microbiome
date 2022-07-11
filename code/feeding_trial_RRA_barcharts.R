##RRA barchart plotting for woodrat lab geeding experiment caecum 16S data

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(RColorBrewer)
library(phyloseq)     # filtering and utilities for such objects

#load phyloseq data
load("physeq_C.rdata")

#call physeq1 for ease of downstream naming scheme
physeq1 <- physeq_C



#remove potential contaminants
physeq1<- subset_taxa(physeq1, Family != "Mitochondria" | Class != "Chloroplast")

#make group species_diet 
physeq1@sam_data$sp_diet <- paste(physeq1@sam_data$Species, physeq1@sam_data$Diet_treatment, sep="_")

#organize data and plot RRA barcharts

########To plot the desired taxanomic level, change each level between this line and the next string of #s below

y1 <- tax_glom(physeq1, taxrank = 'Phylum') # agglomerate taxa at desired level
y2 = suppressWarnings(merge_samples(y1, "sp_diet")) # merge samples on sample variable of interest; suppress warnings is to prevent from failing
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- suppressWarnings(psmelt(y3)) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Phylum[y4$Abundance < 0.05] <- "Phyla < 5% abund." #rename genera with < .25% abundance

#split back out the location and staph into unique columns
y4 <- cbind(y4, read.table(text=y4$Sample, sep="_", header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
y4$Phylum <- factor(y4$Phylum, levels=rev(unique(y4$Phylum)))

#change order so lepida on left
y4$col1 <- factor(y4$col1, levels = c("N. lepida", "N. bryanti"))


#set color palette to accommodate the number of genera
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(y4$Phylum)))


#plot
p <- ggplot(data=y4, aes(x=col2, y=Abundance, fill=Phylum, Genus =Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill=guide_legend(nrow=5)) + 
  scale_x_discrete(limits=rev(levels(as.factor(y4$col2)))) + facet_wrap(~col1) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  ylab("16S relative read abundance") + xlab("")
p


#Genus
ggsave(plot=p, "../Lab-diet-trial-16S-analysis/figures/phyla_RRA.jpg", width = 14, height =8 , device='jpg', dpi=500)
########


#This is just a quick plotting of data to dummy check my above resulting figures,
#change names and taxonomic levels as needed

clostridium <- subset_taxa(physeq1, Genus=="Allobaculum")

plot_bar(clostridium) + facet_wrap(~sp_diet)
