#Deseq2 analysis for caecum microbiome 

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

library(phyloseq)
library(ggplot2)
library(DESeq2)
library(data.table)


#load in phyloseq object for Caecum data
physeq_C <- readRDS("physeq_C.rds")

#rename to WR_species so it doesn't conflict below
names(physeq_C@sam_data)[names(physeq_C@sam_data) == 'Species'] <- 'WR_species'

sum(sample_sums(physeq_C)) #total number of reads 202388 across 1444 otus


##Rarefy?
#set.seed(666)
#physeq <- rarefy_even_depth(physeq,sample.size = min(sample_sums(physeq)))


physeq <- physeq_C
#make factor levels for diet
physeq@sam_data$Diet_treatment <- as.factor(physeq@sam_data$Diet_treatment)


#make a taxonomic rank that includes both family and genus - this will help when visualizing later
Fam_Genus_species <- paste(tax_table(physeq)[ ,"Class"], tax_table(physeq)[ ,"Family"],tax_table(physeq)[ ,"Genus"], tax_table(physeq)[,"Species"], sep = "_")
tax_table(physeq) <- cbind(tax_table(physeq), Fam_Genus_species)


### Start here for N. lepida between diets ###

#subset to only N. lepida to look for diet related DE taxa
physeq_lep <- subset_samples(physeq, physeq@sam_data$WR_species == "N. lepida")

#set FRCA as reference level for diet treatment such that negative values will be greater in FRCA 
#and positive is greater in PRFA
physeq_lep@sam_data$Diet_treatment = relevel(physeq_lep@sam_data$Diet_treatment, "FRCA")

#create the deseq object and run the function  
diff_abund = phyloseq_to_deseq2(physeq_lep, ~Diet_treatment)
diff_abund = DESeq(diff_abund, test="Wald", fitType ="parametric", sfType = "poscounts") #the poscounts sfType flag is used as we have many 0s in the data

#check the reference levels
resultsNames(diff_abund)


#set contrast to compare each species on their native diet
res = results(diff_abund, contrast = list(c( "Diet_treatment_PRFA_vs_FRCA")))
sigtab = res[which(res$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

#Set up plotting framework

#ggplot settings
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Fam_Genus_species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Fam_Genus_species = factor(as.character(sigtab$Fam_Genus_species), levels=names(x))
#20 total DE OTUs 

#then plot results
diff_abund_plot <- ggplot(sigtab, aes(x=Fam_Genus_species, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=1.0) + 
  annotate("text", x = 8, y=-10, label = "FRCA",size = unit(14, "pt")) +
  annotate("text", x = 8, y=10, label = "PRFA",size = unit(14, "pt")) +
  theme(axis.text.x = element_text(size = 20)) + xlab("Taxa") +
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24))

#positive values are taxa that were more abundant in N. lepida individuals on native PRFA diet, 
#and negative values were more abundant in N. lepida on non-native FRCA diet treatment
print(diff_abund_plot)

ggsave(plot=diff_abund_plot, "../Lab-diet-trial-16S-analysis/figures/lep_deseq_PRFAvFRCA.jpg", width = 16, height =10 , device='jpg', dpi=500)


unit <- subset_taxa(physeq_lep, Fam_Genus_species == "Clostridiaceae_Clostridium_NA")

plot_bar(unit) + facet_wrap(~ Diet_treatment)



### Start here for N. bryanti between diets ###

#subset to only N. lepida to look for diet related DE taxa
physeq_bry <- subset_samples(physeq, physeq@sam_data$WR_species == "N. bryanti")

#set PRFA as reference level for diet treatment such that negative values will be greater in PRFA 
#and positive is greater in FRCA
physeq_bry@sam_data$Diet_treatment = relevel(physeq_bry@sam_data$Diet_treatment, "PRFA")

#create the deseq object and run the function  
diff_abund_bry = phyloseq_to_deseq2(physeq_bry, ~Diet_treatment)
diff_abund_bry = DESeq(diff_abund_bry, test="Wald", fitType ="parametric", sfType = "poscounts") #the poscounts sfType flag is used as we have many 0s in the data

#check the reference levels
resultsNames(diff_abund_bry)


#set contrast to compare each species on their native diet
res_bry = results(diff_abund_bry, contrast = list(c("Diet_treatment_FRCA_vs_PRFA")))
sigtab_bry = res_bry[which(res_bry$padj < 0.05), ]
sigtab_bry = cbind(as(sigtab_bry, "data.frame"), as(tax_table(physeq_bry)[rownames(sigtab_bry), ], "matrix"))

#Set up plotting framework

#ggplot settings
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x_bry = tapply(sigtab_bry$log2FoldChange, sigtab_bry$Phylum, function(x_bry) max(x_bry))
x_bry = sort(x_bry, TRUE)
sigtab_bry$Phylum = factor(as.character(sigtab_bry$Phylum), levels=names(x_bry))
# Genus order
x_bry = tapply(sigtab_bry$log2FoldChange, sigtab_bry$Fam_Genus_species, function(x_bry) max(x_bry))
x_bry = sort(x_bry, TRUE)
sigtab_bry$Fam_Genus_species = factor(as.character(sigtab_bry$Fam_Genus_species), levels=names(x_bry))


#then plot results
diff_abund_bry_plot <- ggplot(sigtab_bry, aes(x=Fam_Genus_species, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() +
  theme(axis.text.x = element_text(size = 20)) + xlab("Taxa") +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=1.0) + 
  annotate("text", x = 1, y=10, label = "FRCA",size = unit(14, "pt")) +
  #annotate("text", x = 1, y=-10, label = "PRFA",size = unit(14, "pt")) +
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24))
#positive values are taxa that were more abundant in N. bryanti individuals on native FRCA diet, 
#and negative values were more abundant in N. bryanti on non-native PRFA diet treatment
print(diff_abund_bry_plot)

ggsave(plot=diff_abund_bry_plot, "../Lab-diet-trial-16S-analysis/figures/bry_deseq_FRCAvPRFA.jpg", width = 14, height =8 , device='jpg', dpi=500)


#idiot check the results

physeq_bry <- transform_sample_counts(physeq_bry, function(x) x/sum(x))

unit <- subset_taxa(physeq, Family == "S24-7")

plot_bar(unit) + facet_wrap(~ WR_species +Diet_treatment)

#save diff expressed microbes to data for WGCNA analysis

sigtab$otu <- rownames(sigtab) #add the otu name for the list below, we will match this for WGCNA

write.table(sigtab$otu, "DESeq_lep_diet.tsv", sep="\t")



