#Boxplots for DE genes for Marjorie

rm(list=ls())

setwd("~/Desktop/lab_trials_16S/DESeq2_Neotoma_transcriptomics/Data")

dir.create("../gene_boxplots")

library(ggplot2)

#read in data
all_counts <- readRDS("all_gene_counts.RDS")
sampleTable <- readRDS("sampleTable.RDS")
sampleTable$Tissue <- factor(sampleTable$Tissue, levels=c("FG", "SI", "C", "L"))


#caecum DEGs
caecum.int.DEGs <- read.table("caecum_DEGs/interaction_result.txt", header=T, row.names = 1)
liver_spec_DEGs <- read.table("liver_DEGs/species_res.txt", header=T, row.names = 1)

#to look for gene families

View(all.counts[grep("NBRY_CYP", row.names(all.counts)),]) #replace SULT with desired gene
View(liver_spec_DEGs[grep("NBRY_CYP\\b[1-3][ab]\\b", row.names(liver_spec_DEGs)),])

#create gene specific objects
#cyps
all_cyps <- t(all_counts[grep("NBRY_CYP", row.names(all_counts)),])
all_cyps <- cbind(all_cyps,sampleTable) #add sample info
all_cyps_w <- all_cyps[,c(119:127,1:118)]

#sults
all_sults <- t(all_counts[grep("NBRY_SULT", row.names(all_counts)),])
all_sults <- cbind(all_sults,sampleTable) #add sample info
all_sults_w <- all_sults[,c(29:37,1:28)]

#GSTs
all_GSTs <- t(all_counts[grep("NBRY_GST", row.names(all_counts)),])
all_GSTs <- cbind(all_GSTs,sampleTable) #add sample info
all_GSTs_w <- all_GSTs[,c(32:40,1:31)]

#ABCs
all_ABCs <- t(all_counts[grep("NBRY_ABC", row.names(all_counts)),])
all_ABCs <- cbind(all_ABCs,sampleTable) #add sample info
all_ABCs_w <- all_ABCs[,c(67:75,1:66)]

#TSTs
all_TSTs <- t(all_counts[grep("NBRY_TST", row.names(all_counts)),])
all_TSTs <- cbind(all_TSTs,sampleTable) #add sample info
all_TSTs_w <- all_TSTs[,c(5:13,1:4)]

#FMOs
all_FMOs <- t(all_counts[grep("NBRY_FMO[0-9]", row.names(all_counts)),])
all_FMOs <- cbind(all_FMOs,sampleTable) #add sample info
all_FMOs_w <- all_FMOs[,c(10:18,1:9)]

#UGTs
all_UGTs <- t(all_counts[grep("NBRY_UGT", row.names(all_counts)),])
all_UGTs <- cbind(all_UGTs,sampleTable) #add sample info
all_UGTs_w <- all_UGTs[,c(28:36,1:27)]



#create a total count value for each gene family 
all_cyps_w$total_cyps <- rowSums(all_cyps_w[,10:127])
all_sults_w$total_sults <- rowSums(all_sults_w[,10:37])
all_GSTs_w$total_GSTs <- rowSums(all_GSTs_w[,10:40])
all_ABCs_w$total_ABCs <- rowSums(all_ABCs_w[,10:75])
all_TSTs_w$total_TSTs <- rowSums(all_TSTs_w[,10:13])
all_FMOs_w$total_FMOs <- rowSums(all_FMOs_w[,10:18])
all_UGTs_w$total_UGTs <- rowSums(all_UGTs_w[,10:36])



#plot these totals for gene families
options(scipen = 10000)
#cyps
cyp_boxplot <- ggplot(data = all_cyps_w, aes(x = Tissue, y=total_cyps, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All Cyps") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#sults
sult_boxplot <- ggplot(data = all_sults_w, aes(x = Tissue, y=total_sults, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All Sults") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#GSTs
GST_boxplot <- ggplot(data = all_GSTs_w, aes(x = Tissue, y=total_GSTs, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All GSTs") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#ABCs
ABC_boxplot <- ggplot(data = all_ABCs_w, aes(x = Tissue, y=total_ABCs, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All ABCs") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#TSTs
TST_boxplot <- ggplot(data = all_TSTs_w, aes(x = Tissue, y=total_TSTs, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All TSTs") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#FMOs
FMO_boxplot <- ggplot(data = all_FMOs_w, aes(x = Tissue, y=total_FMOs, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All FMOs") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)

#UGTs
UGT_boxplot <- ggplot(data = all_UGTs_w, aes(x = Tissue, y=total_UGTs, fill=Species)) +
  geom_boxplot() + theme_bw() + ggtitle("All UGTs") + ylab("Normalized Count") + xlab("Tissue") +
  scale_fill_manual(values= c("forestgreen","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size = 20, face= "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18)) + theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 18)) + theme(axis.text.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size=18))+
  theme(strip.text = element_text(size=18)) +
  facet_wrap(~Diet)


ggsave(plot=UGT_boxplot, "../gene_boxplots/UGT.jpg", width = 14, height =8 , device='jpg', dpi=500)

