##################### In Silico Analysis of Metabolic Effects of Bipolar Disorder on Prefrontal Cortex Identified Altered GABA, Glutamate-Glutamine Cycle, Energy Metabolism and Amino Acid Synthesis Pathways###############

##################### DATA ANALYSIS AND DIFFERENTIAL EXPRESSION ##################

##### Required Packages ######

library(affy)
library(Affyhgu133aExpr)
library(hgu133a.db)
library(limma)
library(ggplot2)
library(stringr)
library(rgl)
library(factoextra)
library(NMF)
library(RColorBrewer)
library(devtools)
library(clusterProfiler)
library(factoextra)
library(expss)
library(org.Hs.eg.db)
library(PCSF)


mydata <- ReadAffy() # CEL FILES ARE DOWNLOADED FROM NCBI GEO (GSE5388)

Annot = data.frame(SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=","), 
                   DESC=sapply(contents(hgu133aGENENAME), paste, collapse=","),
                   ENTREZID=sapply(contents(hgu133aENTREZID), paste, collapse=","),
                   ENSEMBLID=sapply(contents(hgu133aENSEMBL), paste, collapse=",")) ### ANNOTATIONS


data.rma.norm = rma(mydata) ### RMA NORMALIZATION 
rma = exprs(data.rma.norm) 

rownames(rma) <- Annot$SYMBOL ### GENE SYMBOLS AS ROWNAMES

rma <- avereps(rma , ID = rownames(rma))

condition <- c(rep("BP",30) , rep("HLT",31)) ### DEFINITION OF CONDITIONS
colnames(rma) <- c(paste("BP" , c(1:30)),paste("HLT" , c(1:31))) ### CONDITIONS AS COLUMN NAMES OF EXPRESSION MATRIX

pca_all <- prcomp(t(rma)) ### PCA
pca_res <- as.data.frame(pca_all$x[,1:2]) ### FIRST 2 PCA
colnames(pca_res) <- c("PC1" , "PC2")

### DATA MATRIX FOR THE USAGE OF PATIENT INFO DOWNLOADED FROM GSE5388 PAGE (ALSO PROVIDED)
all_data <- read.table("GSE5388_data_matrix.csv" , header = F, sep = "," , row.names = 1)
all_data <- as.matrix(all_data)
condition <- all_data["Condition",]
age <- all_data["Age",]
age <- as.numeric(as.matrix(age))
gender <- all_data["Gender",]
duration <- all_data["Duration",]
duration <- as.numeric(duration)
brain_Side <- all_data["Brain Side",]
li_treat <- all_data["Lithium Treatment",]
Alcohol <- all_data["Alcohol",]
pca_res$cond <- condition
pca_res$gender <- gender
pca_res$side <- brain_Side
pca_res$age <- age
pca_res$duration <- duration
pca_res$lithium <- li_treat
pca_res$alcohol <- Alcohol
pca_res$data <- c(paste("BP" , c(1:30)),paste("HLT" , c(1:31)))

### PCA PLOTS OF PATIENT INFO
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = cond)) + geom_point() + theme_bw() + xlab("PC1 (13.11%)") + ylab("PC2 (10.70%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = gender)) + geom_point() + theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = side)) + geom_point() +  theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = age)) + geom_point() + theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = data)) + geom_point() + theme_bw() + geom_text(aes(label = data)) + theme(legend.position="none") + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = lithium)) + geom_point() + theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = alcohol)) + geom_point() + theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = duration)) + geom_point() + theme_bw() + xlab("PC1 (21.6%)") + ylab("PC2 (11.4%)")


##### FIGURE 1 OF ARTICLE ########

fviz_pca_ind(pca_all, geom.ind = "point", pointshape = 21, 
             pointsize = 3.2, 
             fill.ind = condition, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             col.var = "black",
             repel = TRUE,
             legend.title = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

### SEPERATION OF DATA AS HEALTHY AND BIPOLAR DISORDER

bp_data <- rma[,1:30]
bp_ages <- age[1:30]
hlt_data <- rma[,31:61]
hlt_ages <- age[31:61]

### EXTRACTION OF SELECTED POPULATION
#### SELECTION ARE MADE MANUALLY

selected_bp <- bp_data[,c(28,21,17,24,27,25,14,19,23)]
selected_bp_ages <- bp_ages[c(28,21,17,24,27,25,14,19,23)]
selected_bp_gender <- gender[1:30][c(28,21,17,24,27,25,14,19,23)]
selected_bp_side <- brain_Side[1:30][c(28,21,17,24,27,25,14,19,23)]

selected_hlt <- hlt_data[,c(1,5,6,10,16,18,21,23,28)]
selected_hlt_ages <- hlt_ages[c(1,5,6,10,16,18,21,23,28)]
selected_hlt_gender <- gender[31:61][c(1,5,6,10,16,18,21,23,28)]
selected_hlt_side <- brain_Side[31:61][c(1,5,6,10,16,18,21,23,28)]

exp_data <- cbind(selected_hlt , selected_bp)
ages <- c(selected_hlt_ages , selected_bp_ages)
genders <- c(selected_hlt_gender , selected_bp_gender)
sides <- c(selected_hlt_side , selected_bp_side)

#################### ANALYSIS OF SELECTED POPULATION AND DIFFERENTIAL EXPRESSION ANALYSIS ##########

pca <- prcomp(t(exp_data))
pca_res <- as.data.frame(pca$x[,1:2])
pca_res$cond <-  c(rep("HLT",9) , rep("BP",9))
pca_res$gender <- genders
pca_res$side <- sides
pca_res$age <- ages

ggplot(pca_res , aes(x = PC1 , y = PC2 , color = cond)) + geom_point() + theme_bw()
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = gender)) + geom_point() + theme_bw()
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = side)) + geom_point() + theme_bw()
ggplot(pca_res , aes(x = PC1 , y = PC2 , color = age)) + geom_point() + theme_bw()


############ DIFFERENTIAL EXPRESSION WITH LIMMA ################

design <- cbind(Healthy = c(rep(1, 9), rep(0, 9)), Bipolar = c(rep(0, 9), rep(1, 9)))
fit <- lmFit(exp_data, design)
cont.matrix <- makeContrasts(HealthyvsBipolar = Healthy - Bipolar, levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, number = Inf, coef = "HealthyvsBipolar" , adjust = "BH")

sig_genes <- result[result$adj.P.Val < 0.01,]

############ AVERAGE EXPRESSION OF GENES FOR EACH CONDITION FOR FUTURE USE ############

hlt_avg <- rowMeans(exp_data[,1:9])
bp_avg <- rowMeans(exp_data[,10:18])

avg_data <- cbind(hlt_avg , bp_avg)
rownames(avg_data) <- rownames(exp_data)

############ SAVED FILES WILL BE USED ANALYSES ON GENOME SCALE METABOLIC MODEL IN MATLAB #########

write.table(result , file = "limma_res.txt" , quote = F , row.names = T , col.names = T , sep = "\t" , dec = ",")
write.table(avg_data , file = "avg_exp_bp.txt" , sep = "\t" , row.names = T , col.names = T , quote = FALSE) 

############ ENRICHMENT ANALYSES ###############

###### ALL GENES FROM EXPRESSION DATA USED AS BACKGROUND (GENE UNIVERSE FOR MORE PRECISE RESULTS) #####

gene_universe <- rownames(result)
up_genes <- rownames(sig_genes[sig_genes$logFC > 0,]) # UPREGULATED DE GENES
down_genes <- rownames(sig_genes[sig_genes$logFC < 0,]) # DOWNREGULATED DE GENES

## ENRICHMENT ANALYSIS USING GENEONTOLOGY
go_res_up <- enrichGO(up_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL" , ont = "ALL" , universe = gene_universe)
go_res_down <- enrichGO(down_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL" , ont = "ALL" , universe = gene_universe)

## DOT PLOTS OF ENRICHMENT RESULTS (FIGURE 2 and SUPPLEMENTARY FIGURE 1)

dotplot(go_res_up, showCategory=10, split="ONTOLOGY" , label_format = 100) + facet_grid(ONTOLOGY~. , scale = "free") + theme_gray()
dotplot(go_res_down, showCategory=10, split="ONTOLOGY" , label_format = 100) + facet_grid(ONTOLOGY~. , scale = "free") + theme_gray()



### INTEGRATION WITH PPI USING PCSF ###

FCs <- sig_genes$logFC
names(FCs) <- rownames(sig_genes)

terminals <- FCs
data("STRING")
ppi <- construct_interactome(STRING)
subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)

write_graph(subnet, "ppi_bp.graphml", format = "graphml")
