
#WGCNA analysis for caecum gene and microbiome 

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")


install.packages("RSQLite")
library(RSQLite)
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("WGCNA", force = TRUE)

library("ape")          # read tree file
library("Biostrings")   # read fasta file
library("phyloseq")     # filtering and utilities for such objects
library("biomformat") # perhaps unnecessary 
library(ggplot2)
library(picante) 
library(WGCNA)
library(DESeq2)
library(data.table)
library(tibble)
library(ggforce)

#load in microbe and gene count data

load("DE_gene_microbe_counts.rdata")
gene_microbe_counts

#load the sample metadata
load("df_gene_microbe.rdata")
df

#load physeq data
load("physeq_C.rdata")


#check for excessive missing values and remove outliers
#check for features with too many missing values
gsg = goodSamplesGenes(gene_microbe_counts, verbose = 3)
gsg$allOK #if false, remove below

#function for removal of bad features (from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(gene_microbe_counts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(gene_microbe_counts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  gene_microbe_counts_clean = gene_microbe_counts[gsg$goodSamples, gsg$goodGenes]
}

#note that after removing outlier samples, there may be more featuers that need be removed
#i.e., those that are now all zeros

i <- (colSums(gene_microbe_counts, na.rm=T) != 0)
gene_microbe_counts_clean <- gene_microbiome_counts_clean[, i]

#check for outlier samples and remove

sampleTree = hclust(dist(gene_microbe_counts_clean), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#looks like S023373 may be an outlier, so we can remove

gene_microbe_counts_clean <- gene_microbe_counts[row.names(gene_microbe_counts) != "S023381", , drop = FALSE]
#also remove from our df
df_clean <- df[row.names(df) != "S023381", , drop = FALSE]


###DONE WITH CLEAN UP###

#make a DESeq object to normalize below
dds <- DESeqDataSetFromMatrix (
  countData = t(gene_microbe_counts_clean),
  colData = df_clean,
  design = ~1)

#normalize
gene_microbe_counts_norm <- dds

gene_microbe_counts_norm <- varianceStabilizingTransformation(dds)
norm_counts <- assay(gene_microbe_counts_norm) %>% 
  t()


###set up the WGCNA 

allowWGCNAThreads()

power <- c(c(1:10), seq(from=12, to=30, by=2)) # soft-threshold powers to inspect

sft<- pickSoftThreshold(norm_counts, powerVector = power, verbose = ) #call network topology function

#plot the results and inspect power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="SFT (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=power, col="purple")

#based on this - use power of 5

#construct network and detect modules
network1 <- blockwiseModules(norm_counts, power= 5,
                             TOMType = "unsigned", minModuleSize = 3,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             randomSeed = 666, #hail satan!
                             saveTOMFileBase = "/Users/dannynielsen/Desktop/lab_trials_16S",
                             verbose=3)

#inspect number of modules (top row) detected, and number of nodes (botton row) in each
table(network1$colors)

#save module assignments and eigengene (ME) info for next steps

moduleLabels <- network1$colors
moduleColors <- labels2colors(network1$colors)
MEs <- network1$MEs
genetree <- network1$dendrograms[[1]]
save(moduleLabels, moduleColors, MEs, genetree, file = "/Users/dannynielsen/Desktop/lab_trials_16S/WGCNA_data.RData")

MEsO <- moduleEigengenes(norm_counts, moduleColors)$eigengenes

MEs <- orderMEs(MEsO)


mergedColors <- labels2colors(network1$colors)

plotDendroAndColors(network1$dendrograms[[1]],
                    mergedColors[network1$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang=0.03, 
                    addGuide = TRUE, guideHang = 0.05)

par(cex=1.0)

plotEigengeneNetworks(MEs, "Eigengene heatmap",marHeatmap = c(0,4,1,0),
                      plotDendrograms = FALSE)

plotEigengeneNetworks(MEs, "Eigengene dendrogram", marHeatmap  = c(0,4,1,0),
                      plotHeatmaps = FALSE)


#explore results
module_eigengenes <- network1$MEs

head(module_eigengenes)

#chck out some treatment effects
all.equal(df_clean$sampleID, rownames(module_eigengenes)) #check that rows are same between metadata and modules

#let's make a species_diet grouping variable
df_clean$species_diet <- paste(df_clean$Species, df_clean$Diet_treatment, sep = "_")


#create a design matrix using this new variable
des_mat <- model.matrix(~ df_clean$species_diet)

#run linear models on these, needs transposed version of matrix

fit <- limma::lmFit(t(module_eigengenes), design = des_mat) #run the models
fit <- limma::eBayes(fit) #empircal Bayes to smooth standard errors


stats_df <- limma::topTable(fit, number=ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

module_ME1 <- module_eigengenes %>%
  tibble::rownames_to_column("sampleID") %>%
  
  dplyr::inner_join(df_clean %>%
                      dplyr::select(sampleID, species_diet),
                    by = c("sampleID" = "sampleID"))


#plot these results

mod_plot <- ggplot(
  module_ME1,
  aes(
    x = species_diet,
    y = ME4,
    color = species_diet
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic() + ggtitle ("Module 4")

ggsave(plot=mod_plot, "../Lab-diet-trial-16S-analysis/figures/WGCNA_module4.jpg", width = 8, height =8 , device='jpg', dpi=500)


###What features make up this module?
#We can use this to explore any of the modules, by their number

##add the microbial taxa names to the module
#make tax table a df
tax_df <- data.frame(physeq_C@tax_table, check.names = FALSE)
tax_df$OTU <- rownames(tax_df)
#make a joined phyla, class, family 
tax_df$p_c_f_g <- paste(tax_df$Phylum, tax_df$Class, tax_df$Family, tax_df$Genus, sep ="_")

#make a key that now has this info
gene_module_key <- tibble::enframe(network1$colors, name = "feature", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

#the program automatically puts an "X" in front of any feature ID that started with a number, 
#this is problem for matching microbial feature names, so we need to undo that!
gene_module_key[1:65,]$feature <- gsub("X", "", gene_module_key[1:65,]$feature)


#use match to get our taxa names
gene_module_key$microbial_taxa <- tax_df$p_c_f_g[match(gene_module_key$feature, tax_df$OTU)]


#we can now extract information for any module of interest
module1_features <- gene_module_key %>%
  dplyr::filter(module == "ME1")

module2_features <- gene_module_key %>%
  dplyr::filter(module == "ME2")

module3_features <- gene_module_key %>%
  dplyr::filter(module == "ME3")

module4_features <- gene_module_key %>%
  dplyr::filter(module == "ME4")



