library(tidyverse)
library(factoextra)
library(pheatmap)

#upload the txt files for RPKMs

SGThy <- read.table("GSM1415670_AGSC_Thy.txt", header = TRUE)
SGKit <- read.table("GSM1415671_AGSC_Kit.txt", header = TRUE)

SCrep1 <- read.table("GSM1202727_SC_RNAseq_Rep1.txt", header = TRUE)
SCrep2 <- read.table("GSM1202728_SC_RNAseq_Rep2.txt", header = TRUE)
SCrep3 <- read.table("GSM1202729_SC_RNAseq_Rep3.txt", header = TRUE)
SCrep4 <- read.table("GSM1202730_SC_RNAseq_Rep4.txt", header = TRUE)
SCrep5 <- read.table("GSM1202731_SC_RNAseq_Rep5.txt", header = TRUE)

STrep1 <- read.table("GSM1202732_ST_RNAseq_Rep1.txt", header = TRUE)
STrep2 <- read.table("GSM1202733_ST_RNAseq_Rep2.txt", header = TRUE)
STrep3 <- read.table("GSM1202734_ST_RNAseq_Rep3.txt", header = TRUE)
STrep4 <- read.table("GSM1202735_ST_RNAseq_Rep4.txt", header = TRUE)
STrep5 <- read.table("GSM1202736_ST_RNAseq_Rep5.txt", header = TRUE)

sperm <- read.table("GSM1202737_M_RNAseq sperm.txt", header = TRUE)


#SCThy and SGKit have 28 more genes than SCyte, STid and sperm

dim(SGThy)
dim(SCrep1)

#remove extra genes from SCyte, STid and sperm

#check is the same genes for all tables
"%notin%" <- Negate("%in%")
SCrep1_extra_1 <- which(SCrep1$gene %notin% SGThy$Genes)
SCrep1_extra_2 <- which(SCrep1$gene %notin% SGKit$Genes)

STrep1_extra_1 <- which(STrep1$gene %notin% SGThy$Genes)
STrep1_extra_2 <- which(STrep1$gene %notin% SGKit$Genes)

sperm_extra_1 <- which(sperm$gene %notin% SGThy$Genes)
sperm_extra_2 <- which(sperm$gene %notin% SGKit$Genes)

#checking that SCrep1, STrep1 and sperm extra genes are the same
#you can perform this check to rep2-5 for SC and ST
identical(SCrep1_extra_1, SCrep1_extra_2,
          STrep1_extra_1,STrep1_extra_2,
          sperm_extra_1,sperm_extra_2)


#see extra genes to remove in SC, ST and sperm
extra_genes <- SCrep1$gene[SCrep1_extra_1]
extra_genes

#remove extra genes from SC, ST and sperm
SCrep1 <- filter(SCrep1, !SCrep1$gene %in% extra_genes)
SCrep2 <- filter(SCrep2, !SCrep2$gene %in% extra_genes)
SCrep3 <- filter(SCrep3, !SCrep3$gene %in% extra_genes)
SCrep4 <- filter(SCrep4, !SCrep4$gene %in% extra_genes)
SCrep5 <- filter(SCrep5, !SCrep5$gene %in% extra_genes)

STrep1 <- filter(STrep1, !STrep1$gene %in% extra_genes)
STrep2 <- filter(STrep2, !STrep2$gene %in% extra_genes)
STrep3 <- filter(STrep3, !STrep3$gene %in% extra_genes)
STrep4 <- filter(STrep4, !STrep4$gene %in% extra_genes)
STrep5 <- filter(STrep5, !STrep5$gene %in% extra_genes)

sperm <- filter(sperm, !sperm$gene %in% extra_genes)

# We calculate mean expression for SC and ST samples

SC <- rowMeans(cbind(SCrep1[,2], 
                         SCrep2[,2],
                         SCrep3[,2],
                         SCrep4[,2],
                         SCrep5[,2]))

ST <- rowMeans(cbind(STrep1[,2],
                         STrep2[,2],
                         STrep3[,2],
                         STrep4[,2],
                         STrep5[,2]))

# we make the RPKM matrix for all the cells populations

RPKM <- cbind(SGThy[,2],
                 SGKit[,2],
                 SC,
                 ST,
                 sperm[,2])

colnames(RPKM) <- c("SGThy1","SGKit","SCyte","STid","Sperm")
rownames(RPKM) <- SGThy$Genes

head(RPKM)

#look for gene of interest
gene <- "Scml2"
RPKM[rownames(RPKM)==gene,]

# calculate z scores
RPKM_zsc <- t(apply(RPKM,1,scale))
colnames(RPKM_zsc) <- colnames(RPKM)
head(RPKM_zsc)

#remove NAs
RPKM_zsc <- na.omit(RPKM_zsc)
dim(RPKM_zsc)

#apply kmean algorithm, selecting K=5 clusters
set.seed(1)
kmeans_all <- kmeans(RPKM_zsc,5, nstart = 25)

#add cluster info to RPKM tables
RPKM_zsc <- cbind(RPKM_zsc, cluster = kmeans_all$cluster)
head(RPKM_zsc)

#verify
kmeans_all$centers
kmeans_all$size

#graph clusters in a PCA plot
fviz_cluster(kmeans_all, 
             data = RPKM_zsc,
             ellipse.type = "euclid", # concentration elipse
             star.plot = FALSE, 
             repel = TRUE,
             label.select=c("Acr","Tnp2","Prm2","Acrbp","Sycp1","Sycp2","Sycp3","Piwil1","Tnp1","Prm1","Pgk2","Izumo1","Gfra1",
                            "Fgfr3","Nanos2","Nanos3","Thy1","Dazl","Kit","Stra8","H2afx","Spa17","Mlh1",
                            "Dkkl1","Top2a","Cd9","Ddx4","Aldh2", "Acrv1","Fgfr3","H2afv","Hist1h2ba"),
             pointsize = 0,
             ggtheme = theme_minimal())

#check cluster with marker genes
RPKM_zsc[which(rownames(RPKM_zsc)=="Thy1"),] #cluster 3 is SGThy1
RPKM_zsc[which(rownames(RPKM_zsc)=="Kit"),] #cluster 1 is SGKit
RPKM_zsc[which(rownames(RPKM_zsc)=="Sycp1"),] #cluster 5 is SCyte
RPKM_zsc[which(rownames(RPKM_zsc)=="Tnp1"),] # cluster 2 is STid
RPKM_zsc[which(rownames(RPKM_zsc)=="Prm2"),] # cluster 4 is sperm

#create a vector to renumber clusters so they are in order fron SGThy1 to sperm
cell_cluster <- as.numeric(ifelse(RPKM_zsc[,6]=="1", "2", 
                       ifelse(RPKM_zsc[,6]=="2", "4",
                              ifelse(RPKM_zsc[,6]=="3", "1",
                                     ifelse(RPKM_zsc[,6]=="4","5","3")))))

RPKM_zsc <- cbind(RPKM_zsc,cell_cluster) # add cell cluster column
RPKM_zsc <- RPKM_zsc[,-6] # remove original cluster numeration column

head(RPKM_zsc)

# upload list for epigenetic enzymes and marker genes
genelist <- read.csv("epigenetic enzymes and MK genes.csv", header=TRUE)

# create a table with only the genes of interest: enzymes and MarKers
Enz_MK_zsc <- RPKM_zsc[rownames(RPKM_zsc) %in% genelist$Symbol,] #matrix with z scores
Enz_MK_RPKM <- RPKM[rownames(RPKM) %in% genelist$Symbol,] # matrix with RPKM values
Enz_MK_RPKM <- cbind(Enz_MK_RPKM,cell_cluster=Enz_MK_zsc[,6]) #add cluster info to RPKM matrix

#convert to data frame
Enz_MK_zsc <- as.data.frame(Enz_MK_zsc)
Enz_MK_RPKM <- as.data.frame(Enz_MK_RPKM)

#match the order of enzymes in the two tables
Enz_MK_zsc <- Enz_MK_zsc[match(genelist$Symbol,rownames(Enz_MK_zsc)),]
Enz_MK_RPKM <- Enz_MK_RPKM[match(genelist$Symbol,rownames(Enz_MK_RPKM)),]

# Bind enzymes ans MK info to z score and RPKM matrices
Enz_MK_zsc <- cbind(Enz_MK_zsc, genelist[,2:12]) 
Enz_MK_RPKM <- cbind(Enz_MK_RPKM, genelist[,2:12])

head(Enz_MK_zsc)
head(Enz_MK_RPKM)

save(Enz_MK_zsc, file="Enzymes MKs zscore matrix.Rdata")
save(Enz_MK_RPKM, file="Enzymes MKs RPKM matrix.Rdata")

#make a heatmap with enzymes and marker genes
figure2_plot <- Enz_MK_zsc[order(Enz_MK_zsc[,"cell_cluster"]),] %>% arrange("Enzyme")
rownames(figure2_plot) <- figure2_plot$Name

pheatmap(figure2_plot[1:5], fontsize = 5, cluster_cols = FALSE, cluster_rows = FALSE)
