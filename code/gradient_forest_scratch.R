rm(list=ls())

setwd("~/Library/CloudStorage/GoogleDrive-dannynielsen@nevada.unr.edu/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#load and install packages
#install phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")

install.packages("extendedForest", repos="http://R-Forge.R-project.org")
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

library(phyloseq)
library(qiime2R)
library(randomForest)
library(extendedForest)
library(gradientForest)
library(vegan)

BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)



#load metadata

all_metadata <- readRDS("lab_trial_final_metadata.RDS")

#read in the DEseq objects for each species alignment
#these will provide the normalized counts for all ~22000 genes
#bryanti aligned counts
bry_C <- readRDS("bry_caecum_dds.RDS")
bry_FG <- readRDS("bry_foregut_dds.RDS")
bry_L <- readRDS("bry_liver_dds.RDS")
bry_SI <- readRDS("bry_SI_dds.RDS")
#create dataset with all counts across tissues for bryanti
bry_all_tissues <- cbind(bry_C@assays@data$counts,bry_FG@assays@data$counts,
                         bry_L@assays@data$counts,bry_SI@assays@data$counts)

#lepida aligned counts
lep_C <- readRDS("lep_C_dds.RDS")
lep_FG <- readRDS("lep_FG_dds.RDS")
lep_L <- readRDS("lep_liver_dds.RDS")
lep_SI <- readRDS("lep_SI_dds.RDS")
#create dataset with all counts across tissues for lepida
lep_all_tissues <- cbind(lep_C@assays@data$counts,lep_FG@assays@data$counts,
                         lep_L@assays@data$counts,lep_SI@assays@data$counts)


#load the lep and bry C interaction DEG genes
bry_C_int_DEGs <- read.csv("bry_aligned_results/bry_C_interaction_DEGs.csv", header=T)
lep_C_int_DEGs <- read.csv("lep_aligned_results/lep_C_interaction_DEGs.csv", header=T)


#load sample metadata

load("df_gene_microbe.rdata")
df <- df[order(rownames(df)),]
df$sp_diet <- paste(df$Species, df$Diet_treatment, sep="_")
df$max_dose <- all_metadata$Max_Dose.[match(df$WR_ID, all_metadata$ID)]
saveRDS(df, "lab_trial_final_metadata.RDS")


#create subset dfs for each sp-diet treatment group
lep_frca_meta <- subset(df, df$Species=="N. lepida" & df$Diet_treatment=="FRCA")

#load the WGCNA module data

network_mods <- read.csv("WGCNA_final_modules.csv", header= TRUE)

#load Picrust KEGG pathways
kegg <- readRDS("KEGG_2_pathways.RDS")
kegg <- t(kegg)

#load in qiime data to phyloseq - has both FG and C counts



#create phyloseq object from qiime 
physeq <- qza_to_phyloseq(features="qiime_output/table.qza",tree="qiime_output/fasttree-tree-rooted.qza",
                          taxonomy = "qiime_output/taxonomy.qza", metadata = "qiime_output/C_FG_metadata.txt")

sum(sample_sums(physeq@otu_table)) #635,508


#make separate for C and FG samples
physeq_C <- subset_samples(physeq, Gut_region=="Caecum")
physeq_FG <- subset_samples(physeq, Gut_region=="Foregut")

#remove singletons
physeq_C <- prune_taxa(taxa_sums(physeq_C) > 1, physeq_C)
physeq_FG <- prune_taxa(taxa_sums(physeq_FG) > 1, physeq_FG)

#get microbe taxa
tax_df <- data.frame(physeq_C@tax_table, check.names = FALSE)
tax_df$OTU <- rownames(tax_df)
#make a joined phyla, class, family 
tax_df$p_f <- paste(tax_df$Phylum, tax_df$Family, sep ="_")

  
#convert to RRA
physeq_C_RRA <- transform_sample_counts(physeq_C, function(x) x/sum(x)) #converts to relative abundance

physeq_C_RRA_genus <- tax_glom(physeq_C_RRA, taxrank = "Genus")


#get microbe taxa
tax_df <- data.frame(physeq_C_RRA_genus@tax_table, check.names = FALSE)
tax_df$OTU <- rownames(tax_df)



microbe_data <- as.data.frame(physeq_C_RRA_genus@otu_table)
microbe_data$otu <- rownames(microbe_data)
microbe_data$genus <- tax_df$Genus[match(rownames(microbe_data),tax_df$OTU)]
microbe_data <- subset(microbe_data, microbe_data$genus != "uncultured")
rownames(microbe_data) <- microbe_data[,31]
microbe_data <- t(microbe_data[,1:29])


#subset the input and microbe datasets 
lep_frca_input <- subset(input_data, rownames(input_data) %in% rownames(lep_frca_meta))
lep_frca_microbe <- subset(microbe_data, rownames(microbe_data) %in% rownames(lep_frca_meta))

kegg_lep_frca <- subset(kegg, rownames(kegg) %in% rownames(lep_frca_meta))

#convert kegg counts to RRA
kegg_lep_frca <- kegg_lep_frca[,1:39]/rowSums(kegg_lep_frca)

#reshape the module DF to wide so that each module value is in a single column
network_mods_wide <- reshape2::dcast(network_mods, sampleID~variable, value.var = "value")
rownames(network_mods_wide) <- network_mods_wide[,1]
network_mods_wide[,1] <- NULL

#subset networks to only lep_frca
network_mods_wide_lep_frca <- subset(network_mods_wide, rownames(network_mods_wide) %in% rownames(lep_frca_meta))

#comibine the two DFs
gene_microbe <- cbind(lep_frca_input,kegg_lep_frca,lep_frca_microbe, network_mods_wide_lep_frca)
gene_microbe_filt <- gene_microbe[,colSums(gene_microbe !=0) > 0] #remove any that sum to 0


#check that rows match up
all.equal(rownames(lep_frca_input), rownames(lep_frca_microbe), rownames(gene_microbe_filt),
          rownames(network_mods_wide_lep_frca))


#make the species and predictor variables for gradient forest 
specs <- colnames(gene_microbe_filt[grep("^NBRY", colnames(gene_microbe_filt))])
preds.picrust<- colnames(gene_microbe_filt[,345:383])
preds.microbes <- colnames(gene_microbe_filt[,384:464])
preds.WGCNAmods <- colnames(gene_microbe_filt[,465:475])

#mantel
gene_dist <- vegdist(lep_frca_input, method = "bray", binary = FALSE)
microbe_dist <- vegdist(lep_frca_microbe, method= "bray", binary=FALSE)
picrust_dist <- vegdist(gene_microbe_filt_norm[grep("^NBRY", colnames(gene_microbe_filt_norm), invert = TRUE)], method="bray", binary=FALSE)
mods_dist <- vegdist(gene_microbe_filt_norm[,465:475], "bray")

mantel.test <- mantel(gene_dist, microbe_dist, method="pearson", permutations =  999)
plot(gene_dist,microbe_dist)

mantel.test.picrust <- mantel(gene_dist, picrust_dist, method="pearson", permutations =  999)
plot(gene_dist,picrust_dist)

mantel.test.modules <- mantel(microbe_dist, mods_dist, method="pearson", permutations = 999)
plot(microbe_dist,mods_dist)

#Gradient Forest
g_forest <- gradientForest(gene_microbe_filt, predictor.vars = preds.picrust, response.vars= specs, ntree = 1000, check.names = FALSE)

  
g_forest$species.pos.rsq #number of "species" (in this case, genes) for which the variables
#(in this case, kegg pathways) have some predictive power


####Play around with KEGG and GO with gene lists
g_forest_genes <- imp.sub$predictors
g_forest_genes <- rownames(data.frame(g_forest$result))
g_forest_genes <- gsub("NBRY_","",g_forest_genes)
g_forest_genes <- gsub("_[0-9]{4}", "", g_forest_genes)

formatted_gene_names <- sapply(g_forest_genes, function(name) {
  formatted <- tolower(substring(name, 2))
  paste0(toupper(substr(name, 1, 1)), formatted)
})

gene_symbols <- formatted_gene_names

keys_to_map_to <- "ENTREZID"
entrez_ids <- select(org.Mm.eg.db, keys = gene_symbols, columns = "ENTREZID", keytype = "SYMBOL")

GO <- enrichGO(gene = as.numeric(entrez_ids$ENTREZID), 'org.Mm.eg.db', pvalueCutoff = 0.05)

#####

g_forest$result
g_forest$res.u


g_forest$overall.imp
g_forest$overall.imp2
g_forest$res
g_forest$imp.rsq #look into this more...



plot(g_forest, plot.type = "Overall.Importance")

barplot(g_forest$overall.imp)


plot(g_forest, plot.type = "Cumulative.Importance")

plot(g_forest, plot.type = "Split.Density")
plot(g_forest, plot.type = "Performance")


#pca

library(DESeq2)

ds <- DESeqDataSetFromMatrix(countData = t(input_data),
                       colData = df,
                       design= ~ sp_diet)

#mantel
gene_dist <- vegdist(input_data, method = "bray")
microbe_dist <- vegdist(microbe_data, method="bray")

CEN = scale(assay(ds), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)

pca$x


biplot <- factoextra::fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 25), #draw top 25 arrows
                          addEllipses = TRUE,
                          habillage = ds$sp_diet,
                          col.ind = ds$sp_diet,
                          ellipse.level=0.95,
                          palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                          geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                          ind.shape = ds$sp_diet,
                          ind.fill = ds$sp_diet,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Microbiome of Woodrat Caecum") +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))


#do Random Forest to see what biomarkers best predict experimental treatment group

df_lep_frca <- subset(df,df$sp_diet=="N. lepida_FRCA")
df_lep_frca$max_dose <- as.factor(df_lep_frca$max_dose)
#create df of predictors with OTUs as columns and samples are rows - these are all for caecum RNA
predictors.microbes <- microbe_data
predictors.genes <- subset(input_data, rownames(input_data) %in% rownames(df_lep_frca))
predictors.genes <- input_data
predictors.modules <- subset(network_mods_wide, rownames(network_mods_wide) %in% rownames(df_lep_frca))
predictors.picrust <- subset(kegg, rownames(kegg) %in% rownames(input_data))

df_lep_frca <- df_lep_frca[order(rownames(df_lep_frca)),] 

#remove columns that sum to 0
predictors.genes <- predictors.genes[,colSums(predictors.genes !=0) > 0]


#create response variables
response <- as.factor(df$sp_diet) #sp_diet treatment group
response_maxdose <- as.factor(df$max_dose)

all.equal(rownames(df_lep_frca), rownames(predictors.genes), rownames(predictors.modules))



#combine to one df
RF.dataframe.genes <- data.frame(response, response_maxdose, predictors.genes) #original
RF.dataframe.genes.sub <- RF.dataframe.genes[,!colnames(RF.dataframe.genes) %in% rownames(imp.sub)]
RF.dataframe.microbes <- data.frame(response, predictors.microbes) #original
RF.dataframe.modules <- data.frame(response_maxdose, predictors.modules) #original
RF.dataframe.picrust <- data.frame(response, predictors.picrust) #original

#subset to lepida only
RF.dataframe.genes.lep <- subset(RF.dataframe.genes, response=="N. lepida_FRCA")
RF.dataframe.genes.lep <- RF.dataframe.genes.lep[,2:369]
RF.dataframe.genes.lep <- RF.dataframe.genes.lep[,colSums(RF.dataframe.genes.lep !=0) > 0]


RF.dataframe.microbes.lep <- subset(RF.dataframe.microbes, response=="N. lepida")
RF.dataframe.microbes.lep <- RF.dataframe.microbes.lep[,2:110]

#use randomForest to train and test model using out of bag OOB error rate

set.seed(666)

WW.classify.genes <- randomForest(response_maxdose~., data=RF.dataframe.genes.lep, ntree=10000, proximity=TRUE, importance=TRUE) 

#x=RF.dataframe.genes[,colnames(RF.dataframe.genes)!="response_maxdose"],
                                  #y=RF.dataframe.genes$response_maxdose, ntree=1000, proximity=TRUE, importance=TRUE) #this runs the randomforest model
#overall model OOB error rate - 6.9% 
#Class error - 
            #N. bryanti_FRCA N. bryanti_PRFA N. lepida_FRCA N. lepida_PRFA class.error
#N. bryanti_FRCA        6               1              0              0   0.1428571
#N. bryanti_PRFA        0               7              0              1   0.1250000
#N. lepida_FRCA         0               0              7              0   0.0000000
#N. lepida_PRFA         0               0              0              7   0.0000000

varImpPlot(WW.classify.genes)
partialPlot(WW.classify.genes, RF.dataframe.genes.lep, "NBRY_EPHA2_0001")

#plots
imp <- importance(WW.classify.genes)
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- dplyr::arrange(imp, desc(MeanDecreaseAccuracy))
imp.sort$predictors <- factor(imp.sort$predictors, levels=imp.sort$predictors)
imp.sub <- imp.sort[1:50,]
imp.sub <- subset(imp.sort, imp.sort$MeanDecreaseAccuracy > 4)

ggplot(imp.sub, aes(x=predictors, y=MeanDecreaseAccuracy)) + 
  geom_bar(stat="identity", fill= "indianred") +
  coord_flip()

set.seed(666)

WW.classify.microbes <- randomForest(response~., data=RF.dataframe.microbes, ntree=10000, proximity=TRUE) #this runs the randomforest model
#overall model OOB error rate - 31.03% 
#N. bryanti_FRCA N. bryanti_PRFA N. lepida_FRCA N. lepida_PRFA class.error
#N. bryanti_FRCA        5               2              0              0   0.2857143
#N. bryanti_PRFA        1               4              1              2   0.5000000
#N. lepida_FRCA         0               0              5              2   0.2857143
#N. lepida_PRFA         0               1              0              6   0.1428571

#plots
imp <- importance(WW.classify.microbes)
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- dplyr::arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels=imp.sort$predictors)
imp.10 <- imp.sort[1:10,]

ggplot(imp.10, aes(x=predictors, y=MeanDecreaseGini)) + 
  geom_bar(stat="identity", fill= "indianred") +
  coord_flip()

set.seed(666)

WW.classify.modules <- randomForest(response_maxdose~., data=RF.dataframe.modules, ntree=10000, proximity=TRUE, importance=TRUE) #this runs the randomforest model
#overall model OOB error rate - % 3.45

 #                   N. bryanti_FRCA N. bryanti_PRFA N. lepida_FRCA N. lepida_PRFA class.error
#N. bryanti_FRCA        6               1              0              0   0.1428571
#N. bryanti_PRFA        0               8              0              0   0.0000000
#N. lepida_FRCA         0               0              7              0   0.0000000
#N. lepida_PRFA         0               0              0              7   0.0000000


set.seed(666)

WW.classify.picrust <- randomForest(response~., data=RF.dataframe.picrust, ntree=10000, proximity=TRUE) #this runs the randomforest model
#overall model OOB error rate - % 58.62

#             Confusion matrix:
#N. bryanti_FRCA N. bryanti_PRFA N. lepida_FRCA N. lepida_PRFA class.error
#N. bryanti_FRCA               3               3              1              0   0.5714286
#N. bryanti_PRFA               1               3              2              2   0.6250000
#N. lepida_FRCA                0               4              1              2   0.8571429
#N. lepida_PRFA                1               1              0              5   0.2857143


#look at how Foregut genes and microbes may predict treatment groups

all_genes <- data.frame(t(readRDS("all_gene_counts.RDS")))
#subset this to only include FG samples
all_genes_FG <- all_genes[grep("_C_", rownames(all_genes)),]
rownames(all_genes_FG) <- sub("_C_S[0-9]{1-2}", "", rownames(all_genes_FG))
#need remove ind. that did not have 16S data
all_genes_FG <- all_genes_FG[1:29,]
all_genes_FG <- all_genes_FG[,colSums(all_genes_FG !=0) > 0] #remove any that sum to 0
all_genes_FG <- all_genes_FG[,1:18850]/rowSums(all_genes_FG) #convert to RRA


#generate metadata 
FG_df <- physeq_FG@sam_data
FG_df <- FG_df[order(FG_df$WR_ID),] #sort so that the rows match the gene df above
FG_df$sp_diet <- paste(FG_df$Species, FG_df$diet_treatment, sep="_")
FG_df$max_dose <- all_metadata$Max_Dose.[match(FG_df$WR_ID, all_metadata$ID)]


#confirm matching rows
all.equal(rownames(all_genes_FG), FG_df$WR_ID)

#need make phyloseq df at genus level

#before converting to RRA, remove unwanted taxa from FG phyloseq
physeq_FG <- subset_taxa(physeq_FG, Family != "Mitochondria" & Order != "Chloroplast" 
                         & Kingdom != "d__Eukaryota" & Kingdom != "Unassigned" & Genus != "uncultured")

#get microbe taxa
FG_tax_df <- data.frame(physeq_FG@tax_table, check.names = FALSE)
FG_tax_df$OTU <- rownames(FG_tax_df)
#make a joined phyla, class, family 
FG_tax_df$p_g <- paste(FG_tax_df$Phylum, FG_tax_df$Genus, sep ="_")


#convert to RRA
physeq_FG_RRA <- transform_sample_counts(physeq_FG, function(x) x/sum(x)) #converts to relative abundance

physeq_FG_RRA_genus <- tax_glom(physeq_FG_RRA, taxrank = "Genus")

#get microbe taxa
FG_tax_df <- data.frame(physeq_FG_RRA_genus@tax_table, check.names = FALSE)
FG_tax_df$OTU <- rownames(FG_tax_df)


FG_microbe_data <- as.data.frame(physeq_FG_RRA_genus@otu_table)
FG_microbe_data$otu <- rownames(FG_microbe_data)
FG_microbe_data$genus <- FG_tax_df$Genus[match(rownames(FG_microbe_data),FG_tax_df$OTU)]
rownames(FG_microbe_data) <- FG_microbe_data[,30]
FG_microbe_data <- t(FG_microbe_data[,1:28])

#reorder so matches other DFs
FG_microbe_data <- FG_microbe_data[order(match(rownames(FG_microbe_data), rownames(FG_df))), ]

#check for rownames the same
all.equal(rownames(FG_df), rownames(FG_microbe_data))


#create df of predictors with OTUs as columns and samples are rows - these are all for caecum RNA
FG_predictors.microbes <- FG_microbe_data
FG_predictors.genes <- all_genes_FG


#create response variables
FG_response <- as.factor(FG_df$sp_diet) #sp_diet treatment group
FG_response_maxdose <- as.factor(FG_df$max_dose)


#combine to one df
FG_RF.dataframe.genes <- data.frame(FG_response, FG_predictors.genes) #original
FG_RF.dataframe.microbes <- data.frame(FG_response, FG_predictors.microbes) #original

#do the RF analysis now
#for genes
FG_RF.genes <- randomForest(x=FG_RF.dataframe.genes[,colnames(FG_RF.dataframe.genes)!="FG_response"],
                            y=FG_RF.dataframe.genes$FG_response, ntree=1000, proximity=TRUE, importance=TRUE) #this runs the randomforest model
#OOB estimate of  error rate: 60.71%

#microbes
FG_RF.microbes <- randomForest(FG_response~., data=FG_RF.dataframe.microbes, ntree=1000, proximity=TRUE, importance=TRUE)

varImpPlot(FG_RF.microbes)
varImpPlot(FG_RF.genes)

#PCA plot for microbiome of FG
#construct the deseq object
physeq_FG@sam_data$sp_diet <- paste(physeq_FG@sam_data$Species, 
                                              physeq_FG@sam_data$diet_treatment, sep="_")
physeq_FG@otu_table <- physeq_FG@otu_table + 1 
dds = phyloseq_to_deseq2(physeq_FG, ~0 + sp_diet)

vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)



# perform a PCA on the data in assay(x) for the selected genes
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)



biplot <- fviz_pca_ind(pca, repel = TRUE, var.axes = TRUE, axes = c(1,2), #use fviz_pca_ind to plot without the vectors
                          addEllipses = TRUE,
                          #select.var = list(contrib =2), #draw top 25 arrows
                          habillage = dds$sp_diet,
                          col.ind = dds$Species,
                          ellipse.level=0.95,
                          palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                          geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                          ind.shape = dds$Diet,
                          ind.fill = dds$Species,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Microbiome of Woodrat Foregut") +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))






