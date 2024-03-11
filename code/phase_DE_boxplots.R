

setwd("~/Library/CloudStorage/GoogleDrive-dannynielsen@nevada.unr.edu/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#the list of speciesxdiet DE genes in the caecum
caecum.int.DEGs <- bry_C_int_DEGs

#this is normalized gene expression data for all detox genes in the caecum
C_RF.dataframe.genes.detox.sub

#subset the detox gene data to only those that are in the DE gene df

caecum_detox_DEgenes <- C_RF.dataframe.genes.detox.sub[,c(1,(2:218)[colnames(C_RF.dataframe.genes.detox.sub)[2:218] %in% caecum.int.DEGs$X])]

#melt data

melted_df <- melt(caecum_detox_DEgenes)

#add phase info

melted_df$phase <- ifelse(grepl("cyp", melted_df$variable, ignore.case = TRUE), "I",
                          ifelse(grepl("sult", melted_df$variable, ignore.case = TRUE), "II", "III"))











