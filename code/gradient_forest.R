rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")


load("SxD_DE_caecum_genes.RData") #input_data

#convert gene count data to RRA
input_data <- input_data[,1:367]/rowSums(input_data)


#load sample metadata
load("df_gene_microbe.rdata")
df <- df[order(rownames(df)),]
df$sp_diet <- paste(df$Species, df$Diet_treatment, sep="_")

#create subset dfs for each sp-diet treatment group
lep_frca_meta <- subset(df, df$Species=="N. lepida" & df$Diet_treatment=="FRCA")


#load Picrust KEGG pathways
kegg <- readRDS("KEGG_2_pathways.RDS")
kegg <- t(kegg)

#load in qiime data to phyloseq - has both FG and C counts

library(phyloseq)
library(qiime2R)
library(extendedForest)
library(gradientForest)
library(vegan)

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



#comibine the two DFs
gene_microbe <- cbind(lep_frca_input,kegg_lep_frca)
gene_microbe_filt <- gene_microbe[,colSums(gene_microbe !=0) > 0] #remove any that sum to 0


#make the species and predictor variables for gradient forest 
specs <- colnames(gene_microbe_filt[grep("^NBRY", colnames(gene_microbe_filt))])
preds <- colnames(gene_microbe_filt[grep("^NBRY", colnames(gene_microbe_filt), invert = TRUE)])

#mantel
gene_dist <- vegdist(lep_frca_input, method = "jaccard", binary = TRUE)
microbe_dist <- vegdist(lep_frca_microbe, method= "jaccard", binary=TRUE)

mantel.test <- mantel(gene_dist, microbe_dist, method="pearson", permutations =  999)
plot(gene_dist,microbe_dist)

#Gradient Forest
g_forest <- gradientForest(gene_microbe_filt, predictor.vars = preds, response.vars = specs, ntree = 100, check.names = FALSE)

  
g_forest$species.pos.rsq #number of genes that have predictive power
  
g_forest$result
g_forest$res.u


g_forest$overall.imp
g_forest$overall.imp
g_forest$overall.imp2
g_forest$res
g_forest$imp.rsq

  
length(g_forest$result)
  

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
