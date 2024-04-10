

setwd("~/Library/CloudStorage/GoogleDrive-dannynielsen@nevada.unr.edu/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data/data")

#the list of speciesxdiet DE genes in the caecum
  caecum.int.DEGs <- bry_C_int_DEGs
  liver.spec.DEGS <- rownames(species_bry_L)

#subset liver DEGs to include only manually annotated detox genes
  liver.spec.DEGS <- liver.spec.DEGS[liver.spec.DEGS %in% rownames(bry_detox_genes)]

  bry_L <- readRDS("bryFirst_liver_dds.RDS")
  bry_L_counts <- as.data.frame(bry_L@assays@data$counts) #these few steps could be made into loop for each dataset
  bry_L_counts <- t(bry_L_counts[rowSums(bry_L_counts !=0) > 0,]) #this will remove genes with 0 counts across all samples, giving us only genes that were expressed in caecum, and transpose
  rownames(bry_L_counts) <- gsub("_S.*$", "", rownames(bry_L_counts)) #make rownames just the WR ID
  bry_L_counts <- sweep(bry_L_counts, 1, rowSums(bry_L_counts), '/') #to make Relative abundance
  

#subset the bry_L_counts to only those in liver.spec.DEGs
  bry_L_counts <- as.data.frame(bry_L_counts[,(1:19044)[colnames(bry_L_counts)[1:19044] %in% colnames(caecum_detox_DEgenes[,2:19])]])
  
#add the treatment group info
  bry_L_counts$C_response <- caecum_detox_DEgenes$C_response

#this is normalized gene expression data for all detox genes in the caecum
  C_RF.dataframe.genes.detox.sub

#subset the detox gene data to only those that are in the DE gene df

  caecum_detox_DEgenes <- C_RF.dataframe.genes.detox.sub[,c(1,(2:218)[colnames(C_RF.dataframe.genes.detox.sub)[2:218] %in% caecum.int.DEGs$X])]

#melt data
  melted_df <- melt(caecum_detox_DEgenes)
  melted_L_df <- melt(bry_L_counts)

#remove the unwanted portions from gene names 
  melted_df$variable <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "", melted_df$variable)
  melted_L_df$variable <- gsub("NBRY_|gene-|_Nmacr|_Nbry", "", melted_L_df$variable)

#add phase info

  melted_df$phase <- ifelse(grepl("cyp", melted_df$variable, ignore.case = TRUE), "Phase I",
                          ifelse(grepl("sult", melted_df$variable, ignore.case = TRUE), "Phase II", "Phase III"))

  melted_L_df$phase <- ifelse(grepl("cyp", melted_L_df$variable, ignore.case = TRUE), "Phase I",
                            ifelse(grepl("sult", melted_L_df$variable, ignore.case = TRUE), "Phase II", "Phase III"))
  
  
phase_DE_plot <- ggplot(melted_L_df, aes(x=variable, y=value, fill=C_response, alpha=C_response)) + 
  geom_boxplot() +
  facet_wrap(~phase, scales="free") +
  #scale_color_manual(values = c("forestgreen","forestgreen", "maroon", "maroon")) +
  scale_fill_manual(values = c("forestgreen", "forestgreen", "maroon", "maroon")) +
  scale_y_log10() +
  scale_alpha_manual(values=c(1, 0.5, 0.5,1)) +  
  guides(
    fill = guide_legend(title = NULL, reverse = TRUE),
    color = guide_legend(title = NULL, reverse = TRUE),
    alpha = guide_legend(title = NULL, reverse = TRUE)
  ) + theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  coord_flip() + ylab("Log transformed count") + xlab("")
  

ggsave(plot=phase_DE_plot, "../Lab-diet-trial-16S-analysis/figures/phase_DE_liver_plot.jpg", width = 12, height =10 , device='jpg', dpi=500)







