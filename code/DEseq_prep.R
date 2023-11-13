#Code to take the HTSeq counts for each lepida and bryanti and make DESeq objects 
#and other for other downstream analysis

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

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
library(ggplot2)

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
BiocManager::install("tximport")
library(tximport)
library(tximportData)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
BiocManager::install("limma")



#Choose which tissue to analyze
#directory <- "bry_gene_counts/caecum/"
#directory <- "bry_gene_counts/foregut/"
#directory <- "bry_gene_counts/liver/"
#directory <- "bry_gene_counts/smallintestine"
directory <- "salmon_map_counts/cecum"

#files from salmon

salmon_files <- file.path("salmon_map_counts/cecum", sampleTable$salmon_sample$files)
txi.cecum <- tximport(salmon_files, type="salmon", txOut = TRUE, importer=read.delim)


#get  names of htseq-count output files
files <-  list.files(directory, pattern = ".sf")

#sample metadata
metadata <- readRDS("lab_trial_final_metadata.RDS")

#individual tissue metadata
#choose proper match to input samples
sampleTable <- read.table("Meta_data/cecum_meta.txt", header = TRUE)
sampleTable$salmon_sample <- cbind(sampleTable,files)

sampleTable$Diet <- as.factor(sampleTable$Diet)
levels(sampleTable$Diet) <- c("PRFA", "FRCA")

#build group variable to facilitate contrasts (rather than specifying interactions model)
#approach detailed in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
sampleTable$group <- factor(paste0(sampleTable$Species,"_",sampleTable$Diet))


#build dds object for DESeq2, design is a linear model formula given variables to test
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ 0 + group)

#build deseq from tximport
salmond_dds <- DESeqDataSetFromTximport(txi.cecum, sampleTable, ~ 0 + group)



#Run DESeq
dds <- DESeq(dds)
resultsNames(dds)

#save dds as identifiable R object to load in downstream analysis
saveRDS(dds, file="bry_caecum_dds.RDS")


vsd <- vst(salmond_dds, blind=FALSE)
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)

#remove the NBRY_ portion of gene name from the dds counts df
rownames(pca$rotation) <- gsub("NBRY_", "", rownames(pca$rotation))

# visualize
biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 25), #draw top 25 arrows
                          #select.var = list(name = c("NBRY_SULT2B1_0001","NBRY_APOA4_0001", "NBRY_FADS3_0001", "NBRY_SULT2A1_0001")),  #alternative to draw specific substitution loadings
                          addEllipses = TRUE,
                          habillage = salmond_dds$group,
                          col.ind = salmond_dds$group,
                          ellipse.level=0.95,
                          palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                          geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                          ind.shape = salmond_dds$group,
                          ind.fill = salmond_dds$group,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Woodrat Caecum") +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))


