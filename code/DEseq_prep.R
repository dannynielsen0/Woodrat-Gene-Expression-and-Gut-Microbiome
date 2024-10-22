#Code to take the HTSeq, or Salmon, counts for each lepida and bryanti and make DESeq objects 
#and other for other downstream analysis

rm(list=ls())

setwd("/Volumes/MatocqLab/Danny/google_drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

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

BiocManager::install("EnhancedVolcano", force = TRUE)
library(EnhancedVolcano)
BiocManager::install("limma")

remove.packages("ggrepel")
install.packages("ggrepel", dependencies = TRUE)

update.packages(ask = FALSE, dependencies = TRUE)

#save load workspace
save.image(file="DEseq_prep.RData")
load("DEseq_prep.RData")

#From the bryanti first pass alighments and counts, choose which tissue to analyze, change bry to lep for lepida aligned data
directory <- "bryFirst_counts_files/caecum/"
#directory <- "bryFirst_counts_files/foregut/"
#directory <- "bryFirst_counts_files/liver/"
#directory <- "bryFirst_counts_files/small_intestine"
#directory <- "bryFirst_counts_files/all_bry_counts"

#From the lepida second pass alighments and counts, choose which tissue to analyze, change bry to lep for lepida aligned data
#directory <- "lepSecond_counts_files/caecum/"
#directory <- "lepSecond_counts_files/foregut/"
#directory <- "lepSecond_counts_files/liver/"
#directory <- "lepSecond_counts_files/small_intestine"

#get  names of htseq-count output files
files <-  list.files(directory, pattern = ".tsv")

#sample metadata
metadata <- readRDS("lab_trial_final_metadata.RDS")


#individual tissue metadata
#choose proper match to input samples
#sampleTable <- read.table("Meta_data/small_intestine_meta.txt", header = TRUE)
sampleTable <- read.table("Meta_data/cecum_meta.txt", header=TRUE)


#sampleTable$salmon_sample <- cbind(sampleTable,files)

sampleTable$Diet <- as.factor(sampleTable$Diet)
levels(sampleTable$Diet) <- c("PRFA", "FRCA")

#build group variable to facilitate contrasts (rather than specifying interactions model)
#approach detailed in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
sampleTable$group <- factor(paste0(sampleTable$Species,"_",sampleTable$Diet))

#subset the sampleTable to only lepida
sampleTable_lep <- droplevels(subset(sampleTable, sampleTable$Species == "N_lepida"))


#build dds object for DESeq2, design is a linear model formula given variables to test
lep_C <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_lep,
                                  directory = directory,
                                  design= ~ 0 + group)

#build deseq from tximport for salmon-processed data
#salmond_dds <- DESeqDataSetFromTximport(txi.cecum, sampleTable, ~ 0 + group)



#Run DESeq
lep_C_dds <- DESeq(lep_C)

#save dds as identifiable R object to load in downstream analysis
saveRDS(bry_first_foregut_dds, file="bryFirst_fg_dds.RDS")
saveRDS(lep_second_caecum_dds, file="lep_second_caecum_dds.RDS")
saveRDS(bry_first_liver_dds, file="bryFirst_liver_dds.RDS")
saveRDS(bry_first_SI_dds, file="bryFirst_SI_dds.RDS")


#load the deseq objects 
bry_first_liver_dds <- readRDS("bryFirst_liver_dds.RDS")
bry_first_caecum_dds <- readRDS("bryFirst_caecum_dds.RDS")
bry_first_foregut_dds <- readRDS("bryFirst_fg_dds.RDS")
bry_first_SI_dds <- readRDS("bryFirst_SI_dds.RDS")

vsd <- vst(bry_first_liver_dds, blind=FALSE)
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)

#remove the NBRY_ portion of gene name from the dds counts df
  rownames(pca$rotation) <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "", rownames(pca$rotation))

#reorder and give full name to tissue
  bry_first_liver_dds$Tissue <- factor(bry_first_liver_dds$Tissue, levels=c("FG","SI", "C", "L"))
  
  bry_first_liver_dds$group <- factor(bry_first_liver_dds$group, levels=c("N_lepida_PRFA","N_lepida_FRCA","N_bryanti_PRFA","N_bryanti_FRCA"))
  levels(bry_first_liver_dds$group) <- c("Specialist - Home", "Specialist - Away", "Generalist - Away", "Generalist - Home")
  
  levels(bry_first_liver_dds$Tissue) <- c("Foregut","Small Intestine", "Caecum", "Liver")

# visualize a summary pca
biplot <- factoextra::fviz_pca_ind(pca, repel = TRUE, axes = c(1,2),
                          #select.var = list(contrib = 25), #draw top 25 arrows
                          #select.var = list(name = c("Sult2b1-1", "APOA4_0001", "FADS3_0001","Sult2a3-7")),  #alternative to draw specific substitution loadings
                          labelsize=20,
                          addEllipses = TRUE,
                          habillage = bry_first_liver_dds$group,
                          col.ind = bry_first_liver_dds$group,
                          ellipse.level=0.95,
                          palette = c("maroon","pink", "lightgreen", "darkgreen"),
                          geom=c("point"), pointsize = 15,   #change to geom=c("point","text") for sample ID
                          ind.fill = bry_first_liver_dds$group,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Liver") +
  scale_shape_manual(values=c(18,15,17,16)) + 
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 50),
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 30),  # Increase legend text size
        legend.title = element_blank(),
        legend.position = "none",  # Move legend to bottom
        legend.direction = "vertical")
    

ggsave(plot=biplot, "../Lab-diet-trial-16S-analysis/figures/liver_gene_expression_pca.jpg", width = 20, height =15 , device='jpg', dpi=700)






#check results and volcano plot

design <- model.matrix(~0 + group, data = sampleTable_lep)


#compare the average effect of eating RHCA to the average effect of eating PRFA
con.diet <- makeContrasts(PRFAvsFRCA = (groupN_lepida_FRCA - groupN_lepida_PRFA),
                          levels=design)
lep.diet.result <- results(lep_C_dds, contrast = con.diet, alpha = 0.05) #results for caecum, and from bry or lep aligned data

lep.diet.result <- as.data.frame(subset(lep.diet.result, lep.diet.result$padj < 0.05 & abs(lep.diet.result$log2FoldChange) > 2))



#Make volcano plots for each effect
#example below looks at diet effect in N. lepida
# Add a column to define significance for coloring

rownames(lep.diet.result) <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "",rownames(lep.diet.result))

lep.diet.result$Significant <- with(lep.diet.result, 
                                    ifelse(padj < 0.05 & abs(log2FoldChange) > 2, 
                                           ifelse(log2FoldChange > 0, "Upregulated on Preferred Diet", "Upregulated on Non-Preferred Diet"), 
                                           "Not Significant"))


#use detox_pattern from detox_genes.R to set a list of labels to use
#lepida
#use detox_pattern from detox_genes.R to set a list of labels to use
lep.diet.result$label <- ifelse(
  lep.diet.result$padj < 0.05 & 
    abs(lep.diet.result$log2FoldChange) > 2 & 
    grepl(paste(detox_pattern, collapse = "|"), rownames(lep.diet.result), ignore.case = TRUE) & 
    !grepl("pseudo", rownames(lep.diet.result), ignore.case = TRUE),  # Exclude rows with 'pseudo'
  rownames(lep.diet.result), 
  NA
)

# Create new columns for label position to avoid overlap (e.g., move labels slightly up/down from points)
lep.diet.result$label_x <- lep.diet.result$log2FoldChange  # Keep x the same
lep.diet.result$label_y <- lep.diet.result$padj + 10  



#save the number of significant genes
num_lep_sig <- sum(lep.diet.result$Significant != "Not Significant", na.rm = TRUE)


# Create the plot
ggplot(lep.diet.result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant), size = 5) +  # Plot points
  scale_color_manual(values = c("grey","darkgreen", "lightgreen"), na.translate=FALSE) +
  #scale_color_manual(values = c("grey","maroon", "pink"), na.translate=FALSE) +  # Custom colors
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p-value cutoff line
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") + 
  geom_label_repel(
    data = subset(lep.diet.result, !is.na(label)),  # Only for rows with labels
    aes(label = label),  # Use the label column for labeling
    size = 5,  # Adjust label text size
    box.padding = 0.5,  # Adjust padding around the label box
    point.padding = 0.5,  # Adjust padding between points and labels
    segment.color = "black",  # Arrow/line color
    segment.size = 0.5,  # Arrow/line thickness
    segment.linetype = "solid",# Line type (e.g., solid, dashed, dotted)
    min.segment.length = 0,
    max.overlaps = Inf, parse = TRUE # Max number of label overlaps before labels are skipped
  ) +
  # Fold-change cutoff lines
  labs(title = "Caecum",
       x = "log2(Fold Change)",
       y = "-log10(adjusted p-value)") + 
  annotate("text", x = -30, y = 45, label = paste(num_lep_sig, "Differentially \nExpressed Genes"), 
           hjust = 0, vjust = 0, size = 4, color = "black") + xlim(-30,30) + ylim(0,65) +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(plot=last_plot(), "../Lab-diet-trial-16S-analysis/figures/gen_C_genes_subset_volcano.pdf", width = 8, height =8 , device='pdf', dpi=500)






