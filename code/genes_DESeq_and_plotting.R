#This will take the output files from DESeq_prep and do the subsequent deseq analysis and other plotting

rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#load libraries
library(reshape2)
library(factoextra)
library(limma)



#read in the various deseq and other data objects

#bryanti aligned counts
bry_C <- readRDS("bry_caecum_dds.RDS")
bry_FG <- readRDS("bry_foregut_dds.RDS")
bry_L <- readRDS("bry_liver_dds.RDS")
bry_SI <- readRDS("bry_SI_dds.RDS")

#lepida aligned counts
lep_C <- readRDS("lep_C_dds.RDS")
lep_FG <- readRDS("lep_FG_dds.RDS")
lep_L <- readRDS("lep_liver_dds.RDS")
lep_SI <- readRDS("lep_SI_dds.RDS")

#read in the gene names for each
bry_genes <- read.table("bry_gene_names.txt", header = FALSE, sep = "\t")
lep_genes <- read.table("lep_gene_names.txt", header = FALSE, sep = "\t")


##Look at gene (sub)family totals (i.e. cyp2B, Ugt, etc...)
bry_c_all_counts <- bry_C@assays@data$counts
bry_CYP2B <- as.data.frame(t(bry_c_all_counts[grep("_CYP|-cyp", rownames(bry_c_all_counts)),]))
bry_UGT2B <- as.data.frame(t(bry_c_all_counts[grep("_UGT|-Ugt", rownames(bry_c_all_counts)),]))
bry_sult <- as.data.frame(t(bry_c_all_counts[grep("_SULT|-Sult", rownames(bry_c_all_counts)),]))
bry_ABC <- as.data.frame(t(bry_c_all_counts[grep("_ABC|-Abc", rownames(bry_c_all_counts)),]))


lep_c_all_counts <- lep_C@assays@data$counts
lep_CYP2B <- as.data.frame(t(lep_c_all_counts[grep("_CYP|-cyp", rownames(lep_c_all_counts)),]))
lep_UGT2B <- as.data.frame(t(lep_c_all_counts[grep("_UGT|-Ugt", rownames(lep_c_all_counts)),]))
lep_sult <- as.data.frame(t(lep_c_all_counts[grep("_SULT|-Sult", rownames(lep_c_all_counts)),]))
lep_ABC <- as.data.frame(t(lep_c_all_counts[grep("_ABC|-Abc", rownames(lep_c_all_counts)),]))



#plot these out

#create a total cyp count value for bryanti
bry_CYP2B$total_cyp2b <- rowSums(bry_CYP2B)
bry_CYP2B <- cbind(bry_CYP2B,as.data.frame(bry_C@colData@listData)[,1:7])
all_bry_cyp2b_long <- melt(bry_CYP2B)
#UGT
bry_UGT2B$total_ugt2b <- rowSums(bry_UGT2B)
bry_UGT2B <- cbind(bry_UGT2B,as.data.frame(bry_C@colData@listData)[,1:7])
all_bry_ugt2b_long <- melt(bry_UGT2B)
#Sult
bry_sult$total_sult <- rowSums(bry_sult)
bry_sult <- cbind(bry_sult,as.data.frame(bry_C@colData@listData)[,1:7])
#ABC
bry_ABC$total_ABC <- rowSums(bry_ABC)
bry_ABC <- cbind(bry_ABC,as.data.frame(bry_C@colData@listData)[,1:7])


#create a total cyp count value for lepida
lep_CYP2B$total_cyp2b <- rowSums(lep_CYP2B)
lep_CYP2B <- cbind(lep_CYP2B,as.data.frame(lep_C@colData@listData)[,1:7])
all_lep_cyp2b_long <- melt(lep_CYP2B)
#UGT
lep_UGT2B$total_ugt2b <- rowSums(lep_UGT2B)
lep_UGT2B <- cbind(lep_UGT2B,as.data.frame(lep_C@colData@listData)[,1:7])
all_lep_ugt2b_long <- melt(lep_UGT2B)
#Sult
lep_sult$total_sult <- rowSums(lep_sult)
lep_sult <- cbind(lep_sult,as.data.frame(lep_C@colData@listData)[,1:7])
#ABC
lep_ABC$total_ABC <- rowSums(lep_ABC)
lep_ABC <- cbind(lep_ABC,as.data.frame(lep_C@colData@listData)[,1:7])



bry_cyp2b_plot <- ggplot(data = all_bry_cyp2b_long, aes(x = Animal , y=value)) +
  geom_col() + theme_bw() + ggtitle("CYP2Bs - from Bry alignment") + ylab("Normalized Count") + xlab("") + facet_wrap(~variable, scales = "free") +
  # scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))


lep_cyp2b_plot <- ggplot(data = all_lep_cyp2b_long, aes(x = Animal , y=value)) +
  geom_col() + theme_bw() + ggtitle("CYP2Bs - from Lep alignment") + ylab("Normalized Count") + xlab("") + facet_wrap(~variable, scales = "free") +
  # scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

#combine the lep and bry total CYP counts into a df and plot side by side

lep_bry_cyp2b <- cbind(lep_CYP2B[,c(122:128,121)], bry_CYP2B$total_cyp2b)
names(lep_bry_cyp2b)[8] <- "total_Lep_cyp2b"
names(lep_bry_cyp2b)[9] <- "total_Bry_cyp2b"

lep_bry_cyp2b_long <- melt(lep_bry_cyp2b) #melt

cyp2b_compare_plot <- ggplot(data = lep_bry_cyp2b_long , aes(x = Animal , y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') + theme_bw() + ggtitle("All CYPS - from both") + ylab("Normalized Count") + xlab("Individual") + facet_wrap(~group, scales = "free") +
  scale_fill_manual(values= c("maroon","forestgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=cyp2b_compare_plot, "../Lab-diet-trial-16S-analysis/figures/cyp2b_compare.jpg", width = 8, height =6 , device='jpg', dpi=500)


#combine the lep and bry total UGT counts into a df and plot side by side

lep_bry_ugt2b <- cbind(lep_UGT2B[,c(60:66,59)], bry_UGT2B$total_ugt2b)
names(lep_bry_ugt2b)[8] <- "total_Lep_UGT2b"
names(lep_bry_ugt2b)[9] <- "total_Bry_UGT2b"


lep_bry_UGT2b_long <- melt(lep_bry_ugt2b) #melt

UGT2b_compare_plot <- ggplot(data = lep_bry_UGT2b_long , aes(x = Animal , y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') + theme_bw() + ggtitle("All UGTs - from both") + ylab("Normalized Count") + xlab("Individual") + facet_wrap(~group, scales = "free") +
  scale_fill_manual(values= c("maroon","forestgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=UGT2b_compare_plot, "../Lab-diet-trial-16S-analysis/figures/UGT_compare.jpg", width = 8, height =6 , device='jpg', dpi=500)


#combine the lep and bry total Sults counts into a df and plot side by side

lep_bry_sult <- cbind(lep_sult[,c(41:47,40)], bry_sult$total_sult)
names(lep_bry_sult)[8] <- "total_Lep_sult"
names(lep_bry_sult)[9] <- "total_Bry_sult"


lep_bry_sult_long <- melt(lep_bry_sult) #melt

Sult_compare_plot <- ggplot(data = lep_bry_sult_long , aes(x = Animal , y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') + theme_bw() + ggtitle("Sults - from both") + ylab("Normalized Count") + xlab("Individual") + facet_wrap(~group, scales = "free") +
  scale_fill_manual(values= c("maroon","forestgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=Sult_compare_plot, "../Lab-diet-trial-16S-analysis/figures/Sult_compare.jpg", width = 8, height =6 , device='jpg', dpi=500)


#combine the lep and bry total ABCs counts into a df and plot side by side

lep_bry_ABC <- cbind(lep_ABC[,c(78:84,77)], bry_ABC$total_ABC)
names(lep_bry_ABC)[8] <- "total_Lep_sult"
names(lep_bry_ABC)[9] <- "total_Bry_sult"


lep_bry_ABC_long <- melt(lep_bry_ABC) #melt

ABC_compare_plot <- ggplot(data = lep_bry_ABC_long , aes(x = Animal , y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') + theme_bw() + ggtitle("ABCs - from both") + ylab("Normalized Count") + xlab("Individual") + facet_wrap(~group, scales = "free") +
  scale_fill_manual(values= c("maroon","forestgreen")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 18, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 16)) #+
#theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14))

ggsave(plot=ABC_compare_plot, "../Lab-diet-trial-16S-analysis/figures/ABC_compare.jpg", width = 8, height =6 , device='jpg', dpi=500)





#PCA plot of gene expression data

#Plot samples in cluster analysis heatmap
vsd <- vst(bry_C, blind=FALSE)#variance stabilizing transformation
sampleDists <- dist(t(assay(vsd)))


#Plot samples in PCA space with important genes as loading vectors
#Want to use fviz_pca_biplot() instead of default ggplot from DESeq2
# perform a PCA on the data in assay(x) for the selected genes
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)

#remove the NBRY_ portion of gene name from the dds counts df
#rownames(pca$rotation) <- gsub("NBRY_", "", rownames(pca$rotation))

# visualize
biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 25), #draw top 25 arrows
                          #select.var = list(name = c("NBRY_SULT2B1_0001","NBRY_APOA4_0001", "NBRY_FADS3_0001", "NBRY_SULT2A1_0001")),  #alternative to draw specific substitution loadings
                          addEllipses = TRUE,
                          habillage = lep_C$group,
                          col.ind = lep_C$Species,
                          ellipse.level=0.95,
                          palette = c("forestgreen", "darkgreen", "darkred", "maroon"),
                          geom=c("point"), pointsize = 5,   #change to geom=c("point","text") for sample ID
                          ind.shape = lep_C$Diet,
                          ind.fill = lep_C$Species,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Woodrat Caecum") +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))


#Get DESeq results and numbers of DE genes
#first create appropriate design matrix:
sampleTable <- as.data.frame(lep_C@colData@listData)[,1:7] #make the sampleTable object
design <- model.matrix(~0 + group, data = sampleTable)
colnames(design) <- levels(sampleTable$group)


#Set up contrast for interaction effect
#Is the effect of species different depending on diet consumed
con.interaction <- makeContrasts(Interaction = ((N_bryanti_PRFA - N_lepida_PRFA)
                                                - (N_bryanti_FRCA - N_lepida_FRCA)),
                                 levels=design)

#The genes that show a significant diet x species interaction
#i.e. "Is the diet effect different across species?"
interaction_result <- results(lep_C, contrast = con.interaction, alpha = 0.05) #results for caecum, and from bry or lep aligned data
summary(interaction_result) 
int_bry_caecum <- as.data.frame(subset(interaction_result, interaction_result$padj < 0.05 & abs(interaction_result$log2FoldChange) > 2))
#343 DE genes in the caecum on bry aligned data
#277 DE genes in the caecum on lep aligned data
write.table(merge(as.data.frame(subset(interaction_result, interaction_result$padj < 0.05 & abs(interaction_result$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "caecum_DEGs/interaction_result.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)




