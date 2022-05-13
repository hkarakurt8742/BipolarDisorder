#### PLOTTING THE SAMPLING RESULTS ####

## PACKAGES 
 
library(reshape2)
library(ggsci)
library(dplyr)
library(ggplot2)
library(devtools)
library(ggrepel)
library(reshape)
library(tidyverse)
library(expss)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(pals)
library(ActivePathways)

## SAMPLING RESULTS FROM MATLAB

sampling_res <- read.table("sampling_res.txt" , header = T , sep = "\t")

sampling_table <- as.character(sampling_res$Pathway)
sampling_table <- as.data.frame(table(sampling_table))
sampling_table <- sampling_table[order(sampling_table$Freq , decreasing = TRUE),]
colnames(sampling_table) <- c("Pathway" , "Frequency")

## COLOR PALETTE

Pathway <- sampling_table$Pathway
myColors <- polychrome(20)
names(myColors) <- Pathway

## FIGURE 4

ggplot(sampling_table, aes(x = Pathway, y = Frequency, fill = Pathway , colour = Pathway)) +
  geom_bar(stat = "identity", color = "white") +
  geom_label(aes(label = Frequency), color = "black") + theme_bw() + 
  scale_fill_manual(name = "Pathway",values = myColors) +
  theme(axis.title.y=element_blank(),axis.title.x=element_blank()) + guides(fill=FALSE) + theme(text = element_text(size=12),axis.text.x = element_text(angle=90, hjust=1))

## LIST OF SELECTED REACTIONS

selected_rxns <- read.table("selected_rxns.txt" , sep = "\t" , header = T)
selected_rxns_plot_data <- cbind(as.numeric(selected_rxns$Healthy) , as.numeric(selected_rxns$Bipolar.Disorder))
rownames(selected_rxns_plot_data) <- paste(selected_rxns$Formula , selected_rxns$Cell.Type , sep = " - ")
colnames(selected_rxns_plot_data) <- c("Healthy","Bipolar.Disorder")
selected_rxns_plot_data <- as.matrix(selected_rxns_plot_data)
selected_rxns_plot_data[,1] <- round(as.numeric(selected_rxns_plot_data[,1]) , digits = 2)
selected_rxns_plot_data[,2] <- round(as.numeric(selected_rxns_plot_data[,2]) , digits = 2)
selected_rxns_plot_data <- selected_rxns_plot_data[(selected_rxns_plot_data[,1] + selected_rxns_plot_data[,2] > 0),]

df.long2<-reshape2::melt(selected_rxns_plot_data)
colnames(df.long2) <- c("Reaction.Name","Condition","value")
df.long2$value <- round(as.numeric(df.long2$value) , digits = 2)
df.long2 <- df.long2 %>% dplyr::na_if(0)

## FIGURE 5 (WITH ALL RESULTS) ## LAST MODIFICATION WAS DONE MANUALLY FOR VISUAL REASONS

ggplot(df.long2, aes(x = Reaction.Name, y = value, fill = Condition, label = value)) +
  geom_bar(stat = "identity") + xlab("Reaction_ID") + ylab("Reaction Flux") +
  geom_text(size = 2.5, position = position_stack(vjust = 0.5)) + scale_fill_brewer(palette = "Dark2") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5),legend.position="top") +
  theme(axis.text.y = element_text(vjust = 1, hjust=1 , size = 6)) + theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=0.5))



# SAMPLING HYPERGEOMETRIC TEST (TABLE 2)
all_sampling_res <- read.table("all_sampling.txt" , header = T , sep = "\t")
gmt_f <- read.GMT("ims571.gmt")
scores_activepathways <- as.matrix(all_sampling_res$P.Value)
rownames(scores_activepathways) <- all_sampling_res$Reaction
activepathways_results <- ActivePathways(scores_activepathways,gmt_f,correction.method = "BH",cytoscape.file.tag = "sampling_hypergeo")
