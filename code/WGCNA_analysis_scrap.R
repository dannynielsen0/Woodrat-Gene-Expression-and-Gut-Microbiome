
#WGCNA analysis for caecum gene and microbiome 

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")


install.packages("RSQLite")
library(RSQLite)
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("WGCNA", force = TRUE)
install.packages("WGCNA")

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

#load the entire 363 DE genes
load("SxD_DE_caecum_genes.RData")
input_data <- input_data[order(rownames(input_data)),]

#just genes
just_genes <- gene_microbe_counts[,21:65]
#just microbes
just_microbes <- gene_microbe_counts[,1:20]

#load the sample metadata
load("df_gene_microbe.rdata")
df <- df[order(rownames(df)),]

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

sampleTree = hclust(dist(gene_microbe_counts), method = "average");
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

#use the 367 dataset
just_genes <- input_data

#make a DESeq object to normalize below
dds <- DESeqDataSetFromMatrix (
  countData = t(just_genes),
  colData = df,
  design = ~1)

#normalize
norm_counts <- assay(dds) %>% t()
norm_counts[1:29, 1:45] <- as.numeric(norm_counts[1:29, 1:45])

#gene_microbe_counts_norm <- varianceStabilizingTransformation(dds)
#norm_counts <- assay(gene_microbe_counts_norm) %>% 
  #t()


###set up the WGCNA 

allowWGCNAThreads()

power <- c(c(1:10), seq(from=12, to=30, by=2)) # soft-threshold powers to inspect

sft<- pickSoftThreshold(norm_counts, powerVector = power, verbose = TRUE) #call network topology function

#plot the results and inspect power - use 12
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="SFT (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=power, col="purple")

#based on this - use power of 12

#this overrides a cor function from another package that conflicts with the cor function in WGCNA
cor <- WGCNA::cor

#construct network and detect modules
network1 <- WGCNA::blockwiseModules(norm_counts, power= 12,
                             TOMType = "signed", minModuleSize = 5, #signed type finds of positively correlated relationships
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             randomSeed = 666, #hail satan!
                             saveTOMFileBase = "/Users/dannynielsen/Desktop/lab_trials_16S",
                             verbose=3)

#inspect number of modules (top row) detected, and number of nodes (bottom row) in each
table(network1$colors)
summary(network1$colors)


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

plotEigengeneNetworks(MEs, "Eigengene heatmap",marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE)

plotEigengeneNetworks(MEs, "Eigengene dendrogram", marHeatmap  = c(0,4,1,0),
                      plotHeatmaps = FALSE)


#explore results
module_eigengenes <- network1$MEs

head(module_eigengenes)

df_clean <- df

#chck out some treatment effects
all.equal(df_clean$sampleID, rownames(module_eigengenes)) #check that rows are same between metadata and modules

#let's make a species_diet grouping variable
df_clean$species_diet <- paste(df_clean$Species, df_clean$Diet_treatment, sep = "_")


#create a design matrix using this new variable
des_mat <- model.matrix(~ df_clean$species_diet)

#run linear models on these to look at relative expression levels, needs transposed version of matrix
#using methods from:
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html

fit <- limma::lmFit(t(MEs), design = des_mat) #run the models
fit <- limma::eBayes(fit) #empircal Bayes to smooth standard errors


stats_df <- limma::topTable(fit, number=ncol(MEs)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

module_ME1 <- MEs %>%
  tibble::rownames_to_column("sampleID") %>%
  
  dplyr::inner_join(df_clean %>%
                      dplyr::select(sampleID, species_diet),
                    by = c("sampleID" = "sampleID"))

mod_long <- melt(data=module_ME1[,c(13,2:12)])

mod_long <- subset(mod_long, mod_long$variable!="MEgrey")

#plot these results


mod_plot <- ggplot(
  mod_long,
  aes(x=species_diet, y = value, fill=variable)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(outlier.shape = NA) + #coord_flip() +
  # A sina plot to show all of the individual data points
  #ggforce::geom_sina(maxwidth = 0.3) + 
  ylab("Expression Level") + xlab("") + ylim(-0.4,0.4) +
  theme_bw() + scale_fill_manual(values= c("pink", "green", "red", "magenta", "blue",
                                      "brown", "purple", "yellow", "black", "turquoise")) +
  theme(plot.title = element_text(size=22)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) + 
  theme(strip.text.x = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title=element_text(size=20)) +
  guides(fill=guide_legend(title="Module"))


ggsave(plot=mod_plot, "../Lab-diet-trial-16S-analysis/figures/modules_grouped_boxplot.jpg", width = 12, height =8 , device='jpg', dpi=500)



#make a key that now has this info
gene_module_key <- tibble::enframe(network1$colors, name = "feature", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key <- gene_module_key[order(gene_module_key$module),]
write.csv(gene_module_key, "network_modules.csv")


#we can now extract information for any module of interest
module1_features <- gene_module_key %>%
  dplyr::filter(module == "ME1")

module2_features <- gene_module_key %>%
  dplyr::filter(module == "ME2")

module3_features <- gene_module_key %>%
  dplyr::filter(module == "ME3")

module4_features <- gene_module_key %>%
  dplyr::filter(module == "ME4")

module5_features <- gene_module_key %>%
  dplyr::filter(module == "ME5")

module6_features <- gene_module_key %>%
  dplyr::filter(module == "ME6")

module7_features <- gene_module_key %>%
  dplyr::filter(module == "ME7")

module8_features <- gene_module_key %>%
  dplyr::filter(module == "ME8")

module9_features <- gene_module_key %>%
  dplyr::filter(module == "ME9")

module10_features <- gene_module_key %>%
  dplyr::filter(module == "ME10")


ggplot(data=gene_microbe_counts, aes(x=NBRY_ANXA1_0001, y=NBRY_BNIP3_0001)) +
  geom_point() + ggtitle("Same Module")

ggplot(data=gene_microbe_counts, aes(x=NBRY_CYP3A11_0004, y=NBRY_SLC10A2_0001)) +
  geom_point() + ggtitle("Neg. correlated modules")


#get the biological data - i.e., microbe counts
biol <- t(just_microbes)

#subset to just N. lepida on FRCA
df_lep_frca <- subset(df, df$Species =="N. lepida")

biol_lep_frca <- data.frame(t(biol[,1:29]))
biol_lep_frca$samID <- rownames(biol_lep_frca)

#sort order of biol_lep by numeric
biol_lep_frca <- biol_lep_frca[order(biol_lep_frca$samID),]


biol_lep_frca <- subset(biol_lep_frca, biol_lep_frca$samID %in% df_lep_frca$sampleID)
biol_lep_frca <- biol_lep_frca[1:14,]
biol_lep_frca <- t(biol_lep_frca[,1:20])#this should work now for the corr. matrix below


rownames(biol_lep_frca) <- gsub("X", "", rownames(biol_lep_frca))#need to remove the X
biol_lep_frca <- data.frame(unlist(biol_lep_frca))
biol_lep_frca <- biol_lep_frca[,1:14]

#get microbe taxa
tax_df <- data.frame(physeq_C@tax_table, check.names = FALSE)
#make a joined phyla, class, family 
tax_df$p_g <- paste(tax_df$Phylum, tax_df$Genus, sep ="_")

#get taxa
biol_lep_frca$microbial_taxa <- tax_df$p_g[match(rownames(biol_lep_frca), rownames(tax_df))]


#cut MEs down to same lepida only
MEs_lep_frca <- MEs[rownames(MEs) %in% colnames(biol_lep_frca),]
MEs_lep_frca <- MEs_lep_frca[,1:10] #remove the grey 'module', which is the 11th module

#make nGenes and NSamples variables to match the lep only data
nGenes = ncol(just_genes)
nSamples = 14


#module to biological data 
moduleCor_frca = cor(MEs_lep_frca, t(biol_lep_frca[,1:14]), use = "p");

modulePvalue_frca = corPvalueStudent(moduleCor_frca, nSamples);
textMatrix = paste(signif(moduleCor_frca, 2), "\n(",
                   signif(modulePvalue_frca, 1), ")", sep = "");
dim(textMatrix) = dim(moduleCor_frca)

#remove grey module for plot
moduleCor_frca <- subset(moduleCor_frca, rownames(moduleCor_frca)!="MEgrey")

sizeGrWindow(16,12)
par(mar = c(14,14, 4, 4))
labeledHeatmap(Matrix = moduleCor_frca,
               xLabels = biol_lep_frca$microbial_taxa,
               yLabels = names(MEs_lep_frca),
               ySymbols = names(MEs_lep_frca),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,  cex.text = 0.65, zlim = c(-1,1),
               main = expression(paste("Module-gene-microbe relationships in ", italic("N. lepida") , "consuming FRCA")))


##Do same for N. bryanti
#subset to just N. bryanti

df_bry <- subset(df, df$Species =="N. bryanti")

biol_bry <- data.frame(t(biol[,1:29]))
biol_bry$samID <- rownames(biol_bry)

#sort order of biol_lep by numeric
biol_bry <- biol_bry[order(biol_bry$samID),]

biol_bry <- subset(biol_bry, biol_bry$samID %in% df_bry$sampleID)
biol_bry <- biol_bry[1:20,]
biol_bry <- t(biol_bry[,1:20])#this should work now for the corr. matrix below


rownames(biol_bry) <- gsub("X", "", rownames(biol_bry))#need to remove the X
biol_bry <- data.frame(unlist(biol_bry))
biol_bry <- biol_bry[,1:15]

#get microbe taxa
tax_df <- data.frame(physeq_C@tax_table, check.names = FALSE)
tax_df$OTU <- rownames(tax_df)
#make a joined phyla, class, family 
tax_df$p_f <- paste(tax_df$Phylum, tax_df$Family, sep ="_")

#get taxa
biol_bry$microbial_taxa <- tax_df$p_f[match(rownames(biol_bry), tax_df$OTU)]


#cut MEs down to same bryida only
MEs_bry <- MEs[rownames(MEs) %in% colnames(biol_bry),]
MEs_bry <- MEs_bry[,1:10]

#make nGenes and NSamples variables to match the lep only data
nGenes = ncol(just_genes)
nSamples = 15

#module to biological data 
moduleCor = cor(MEs_bry, t(biol_bry[,1:15]), use = "p");

modubryvalue = corPvalueStudent(moduleCor, nSamples);
textMatrix = paste(signif(moduleCor, 2), "\n(",
                   signif(modubryvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleCor)

sizeGrWindow(16,12)
par(mar = c(14,14, 4, 4))
labeledHeatmap(Matrix = moduleCor,
               xLabels = biol_bry$microbial_taxa,
               yLabels = names(MEs_bry),
               ySymbols = names(MEs_bry),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,  cex.text = 0.65, zlim = c(-1,1),
               main = expression(paste("Module-gene-microbe relationships in ", italic("N. bryanti"))))







#Exporting data for network figure
options(stringsAsFactors = TRUE);
#resultied in empty matrix

#save the norm counts from above here
TOM=TOMsimilarityFromExpr(norm_counts, power=12)


softPower <- 12 ;
adjacency <- adjacency(norm_counts, power = softPower) ;
TOM <- TOMsimilarity(adjacency)
modules = c("pink", "green", "red", "magenta", "blue", "brown", "purple",
            "yellow", "black", "turquoise", "grey");
probes = colnames(norm_counts)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Cytogene_edges", ".txt", sep=""),
                               nodeFile = paste("Cytogene_nodes", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modules,
                               nodeAttr = moduleColors[inModule])


###sanity check the co-expression of genes in modules

mod_8 <- module8_features$feature


input_data_check <- data.frame(t(input_data))
input_data_check$gene <- rownames(input_data_check)

input_8 <- subset(input_data_check, input_data_check$gene %in% module8_features$feature)
