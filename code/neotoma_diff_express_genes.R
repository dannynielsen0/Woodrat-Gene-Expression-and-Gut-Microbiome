
rm(list=ls())

setwd("~/Desktop/lab_trials_16S/DESeq2_Neotoma_transcriptomics/Data")

#load R packages for DEG data analysis
library(DESeq2)
library(GenomicFeatures)
library(apeglm)
library(ashr)
library(IHW)
library(PoiClaClu)
BiocManager::install("limma")
library(limma)
library(ape)

#load R packages for data visualization
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(plotly)
library(stats)
library(MASS)
library(reshape2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
BiocManager::install("limma")


#set  working directory and subdirectory with counts
dir.create("caecum_DEGs")
dir.create("FG_DEGs")
dir.create("liver_DEGs")
dir.create("SI_DEGs")

dir.create("caecum_plots")
dir.create("FG_plots")
dir.create("liver_plots")
dir.create("SI_plots")

#Choose which tissue to analyze
directory <- "Counts_files/Cecum/"
#directory <- "Counts_files/Foregut/"
#directory <- "Counts_files/Liver/"
#directory <- "Counts_files/SmallIntestine/"


#get  names of htseq-count output files
files <-  list.files(directory, pattern = ".tsv")

#sample metadata loaded here
#choose proper match to input samples
sampleTable <- read.table("Meta_data/cecum_meta.txt", header = TRUE)
sampleTable$Diet <- as.factor(sampleTable$Diet)
levels(sampleTable$Diet) <- c("PRFA", "FRCA")

#build group variable to facilitate contrasts (rather than specifying interactions model)
#approach detailed in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
sampleTable$group <- factor(paste0(sampleTable$Species,"_",sampleTable$Diet))


#build dds object for DESeq2, design is a linear model formula given variables to test
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ 0 + group)

#Run the differential expression analysis on the dataset
dds <- DESeq(dds)
resultsNames(dds)

#read in table of Name and gene_id for results export
genes <- read.table("gene_names.txt", header = TRUE, sep = "\t",
                    row.names = 1)


#set up contrasts of interest
#first create appropriate design matrix:
design <- model.matrix(~0 + group, data = sampleTable)
colnames(design) <- levels(sampleTable$group)

#compare the average effect of being a bryanti to the average effect of being lepida
con.species<- makeContrasts(lepVSbry = (N_lepida_PRFA + N_lepida_FRCA)/2
                            - (N_bryanti_PRFA + N_bryanti_FRCA)/2,
                            levels=design)

#compare the average effect of eating RHCA to the average effect of eating PRFA
con.diet <- makeContrasts(PRFAvsRHCA = (N_lepida_PRFA + N_bryanti_PRFA)/2
                          - (N_lepida_FRCA + N_bryanti_FRCA)/2,
                          levels=design)

#Is the effect of species different depending on diet consumed
con.interaction <- makeContrasts(Interaction = ((N_bryanti_PRFA - N_lepida_PRFA)
                                                - (N_bryanti_FRCA - N_lepida_FRCA)),
                                 levels=design)


#Get table of average effect of species
#Using FDR cutoff of p < 0.05 and required log2fold change above 2
species_res <- results(dds, contrast=con.species, alpha = 0.05)
summary(species_res)
write.table(merge(as.data.frame(subset(species_res, species_res$padj < 0.05 & abs(species_res$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "liver_DEGs/species_res.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)

#plot most significant gene for species:  
ix = which.min(species_res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group)  )
plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("group"), main=rownames(dds[ix,]), )


#The average effect of diet 
diet_res <- results(dds, contrast = con.diet , alpha = 0.05)
summary(diet_res)
write.table(merge(as.data.frame(subset(diet_res, diet_res$padj < 0.05 & abs(diet_res$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "liver_DEGs/diet_res.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)


#plot most significant gene for diet:  
ix = which.min(diet_res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group)  )
plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("group"), main=rownames(dds[ix,]), )

#figure out the pvalue of the 4th lowest we can then use that to subset and plot 4
diet_res$gene <- rownames(diet_res)
species_res$gene <- rownames(species_res)

int_new <- species_res[order(species_res$padj, decreasing = FALSE), ]  # Order data descending
top4 <- int_new[1:4,]
ix = which(species_res$gene=="NBRY_SULT2A1_0001") # most significant
one <- ix[1]
two <- ix[2]
three <- ix[3]
four <- ix[4]
int_new$gene <- rownames(int_new)
sult <- subset(int_new,int_new$gene=="NBRY_SULT2B1_0001")

barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group))

#then plot these
plot_one <- plotCounts(dds, gene=rownames(dds[one,]), las=2, intgroup=c("group"), main=rownames(dds[one,]), returnData =FALSE)
plot_CYP2D_0007 <- ggplot(data = plot_one, aes(x = group, y=count)) +
  geom_boxplot() + theme_bw() + ggtitle("NBRY_SULT2A1_0001") + ylab("Normalized Count") + xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16)) + theme(axis.title.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=plot_CYP2D_0007, "cecum_plots/SULT2A1.jpg", width = 8, height =6 , device='jpg', dpi=500)
                      
plot_two <- plotCounts(dds, gene=rownames(dds[two,]), las=2, intgroup=c("group"), main=rownames(dds[two,]),returnData = TRUE)
plot_CYP2D_0009 <- ggplot(data = plot_two, aes(x = group, y=count)) +
  geom_boxplot() + theme_bw() + ggtitle("CYP2D_0009") + ylab("Normalized Count") + xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16)) + theme(axis.title.y = element_text(size = 16))+
theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=plot_CYP2D_0009, "CYP2D_0009.jpg", width = 8, height =6 , device='jpg', dpi=500)


plot_three <- plotCounts(dds, gene=rownames(dds[three,]), las=2, intgroup=c("group"), main=rownames(dds[three,]),returnData = TRUE)
plot_CYP2D_0010 <- ggplot(data = plot_three, aes(x = group, y=count)) +
  geom_boxplot() + theme_bw() + ggtitle("CYP2D_0010") + ylab("Normalized Count") + xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16)) + theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=plot_CYP2D_0010, "CYP2D_0010.jpg", width = 8, height =6 , device='jpg', dpi=500)


plot_four <- plotCounts(dds, gene=rownames(dds[four,]), las=2, intgroup=c("group"), main=rownames(dds[four,]),returnData = TRUE )




#The genes that show a significant diet x species interaction
#i.e. "Is the diet effect different across species?"
interaction_result <- results(dds, contrast = con.interaction, alpha = 0.05)
summary(interaction_result)
write.table(merge(as.data.frame(subset(interaction_result, interaction_result$padj < 0.05 & abs(interaction_result$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "liver_DEGs/interaction_result.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)


#plot most significant gene for interaction:  
ix = which.min(interaction_result$padj) # most significant

#figure out the pvalue of the 4th lowest we can then use that to subset and plot 4
int_new <- interaction_result[order(interaction_result$padj, decreasing = FALSE), ]  # Order data descending
top4 <- int_new[1:4,]
ix = which(interaction_result$padj < 0.05 & abs(interaction_result$log2FoldChange) > 2) # most significant
one <- ix[1]
two <- ix[2]
three <- ix[3]
four <- ix[4]

barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group))

#then plot these
plot_one <- plotCounts(dds, gene=rownames(dds[four,]), las=2, intgroup=c("group"), main=rownames(dds[four,]), returnData = TRUE)
plot_DNAH17_001 <- ggplot(data = plot_one, aes(x = group, y=count)) +
  geom_boxplot() + theme_bw() + ggtitle("TET3_0017") + ylab("Normalized Count") + xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 16))

ggsave(plot=plot_DNAH17_001, "TET3_0017.jpg", width = 10, height =6 , device='jpg', dpi=500)


plot_two <- plotCounts(dds, gene=rownames(dds[two,]), las=2, intgroup=c("group"), main=rownames(dds[two,]), )
plot_three <- plotCounts(dds, gene=rownames(dds[three,]), las=2, intgroup=c("group"), main=rownames(dds[three,]), )
plot_four <- plotCounts(dds, gene=rownames(dds[four,]), las=2, intgroup=c("group"), main=rownames(dds[four,]), )



#Plot samples in cluster analysis heatmap
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$Species, vsd$Diet, sep="--")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

ggsave(plot=heatmap, "SI_plots/SI_heatmap.pdf", width = 10, height =10 , device='pdf', dpi=500)



#Plot samples in PCA space with important genes as loading vectors
#Want to use fviz_pca_biplot() instead of default ggplot from DESeq2
# perform a PCA on the data in assay(x) for the selected genes
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)
# visualize
biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                #select.var = list(contrib = 25), #draw top 25 arrows
                select.var = list(name = c("NBRY_SULT2B1_0001","NBRY_APOA4_0001", "NBRY_FADS3_0001", "NBRY_SULT2A1_0001")),  #alternative to draw specific substitution loadings
                addEllipses = TRUE,
                habillage = dds$group,
                col.ind = dds$Species,
                ellipse.level=0.95,
                palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                ind.shape = dds$Diet,
                ind.fill = dds$Species,
                invisible = c( "quali"), #remove enlarged symbol for group mean
                title = "Woodrat Caecum") +
                theme(text = element_text(size = 16),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 16))


ggsave(plot=biplot, "caecum_plots/caecum_biplot_interesting_Vec_bry_sult.pdf", width = 10, height =10 , device='pdf', dpi=500)




#Make volcano plots for each effect
#example below looks at diet effect in N. lepida
diet_volc <- merge(as.data.frame(diet_res), 
                   genes, by="row.names", all.y=FALSE)
volcano <- EnhancedVolcano(diet_volc,
                lab = diet_volc$gene_id,
                title = "Volcano plot",
                subtitle = bquote(italic("Average Diet Effect in N. lepida Small Intestine")),
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-5.5, 5.5),
                #ylim = c(0, -log10(10e-12)),
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 2.0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                drawConnectors = TRUE,
                widthConnectors = 0.25,)


ggsave(plot=volcano, "SI_plots/SI_lep_volcano.pdf", width = 12, height =10 , device='pdf', dpi=500)







