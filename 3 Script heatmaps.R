library(tidyverse)
library(pheatmap)

load("RPKM matrix.Rdata")

#heatmap Cluster1
cluster1 <- Enz_MK %>% filter(cell_cluster=="1")
rownames(cluster1) <-  cluster1$Name

#to sort the enzymes and MK make vector with enzymes in same order that buble plots, and add MK names

#vector with cluster 1 marker genes names
markersC1 <- rownames(cluster1[cluster1$Enzyme=="MK",])
####vector with cluster1 enzymes in the same order as buble plot from cluster 1
#Cluster 1 buble plot names are in SGThy1_cluster1, but this table has duplicated names
#to obtain similar distribution, we remove the first entry of duplicated name
enzymesC1 <-SGThy1_cluster1$Name[-which(duplicated(SGThy1_cluster1$Name,fromLast = TRUE))]

#create vector with enzymes and markers order
order_enzC1 <- c(enzymesC1,markersC1)

#reorder cluster 1 rownames with order vector previously made
cluster1 <- cluster1[match(order_enzC1,rownames(cluster1)),]

annot_rowC1 <- data.frame(Enzyme=cluster1$Enzyme)
rownames(annot_rowC1) <- rownames(cluster1)
annot_rowC1$Enzyme <- factor(annot_rowC1$Enzyme, levels = c("Writer","Eraser","MK")) 

annot_colors <- list(Enzyme=c(Writer="cornflowerblue",Eraser="Red",MK="Chartreuse3"))

pheatmap(cluster1[,1:5], fontsize = 10, 
         cluster_cols = FALSE, cluster_rows = FALSE, 
         cellwidth = 30, cellheight = 15,
         annotation_row = annot_rowC1,
         annotation_colors=annot_colors,
         annotation_names_row = FALSE)

##############################################################################################
##############################################################################################
#heatmap Cluster2
cluster2 <- Enz_MK %>% filter(cell_cluster=="2")
rownames(cluster2) <-  cluster2$Name

#to sort the enzymes and MK make vector with enzymes in same order that buble plots, and add MK names

#vector with cluster 2 marker genes names
markersC2 <- rownames(cluster2[cluster2$Enzyme=="MK",])
####vector with cluster2 enzymes in the same order as buble plot from cluster 2
#Cluster 2 buble plot names are in SGKit_cluster2, but this table has duplicated names
#to obtain similar distribution, we remove the first entry of duplicated name
enzymesC2 <-SGKit_cluster2$Name[-which(duplicated(SGKit_cluster2$Name,fromLast = TRUE))]

#create vector with enzymes and markers order
order_enzC2 <- c(enzymesC2,markersC2)

#reorder cluster 2 rownames with order vector previously made
cluster2 <- cluster2[match(order_enzC2,rownames(cluster2)),]

annot_rowC2 <- data.frame(Enzyme=cluster2$Enzyme)
rownames(annot_rowC2) <- rownames(cluster2)
annot_rowC2$Enzyme <- factor(annot_rowC2$Enzyme, levels = c("Writer","Eraser","MK")) 

annot_colors <- list(Enzyme=c(Writer="cornflowerblue",Eraser="Red",MK="Chartreuse3"))

pheatmap(cluster2[,1:5], fontsize = 10, 
         cluster_cols = FALSE, cluster_rows = FALSE, 
         cellwidth = 30, cellheight = 15,
         annotation_row = annot_rowC2,
         annotation_colors=annot_colors,
         annotation_names_row = FALSE)

##############################################################################################
##############################################################################################
#heatmap Cluster3
cluster3 <- Enz_MK %>% filter(cell_cluster=="3")
rownames(cluster3) <-  cluster3$Name

#to sort the enzymes and MK make vector with enzymes in same order that buble plots, and add MK names

#vector with cluster 3 marker genes names
markersC3 <- rownames(cluster3[cluster3$Enzyme=="MK",])
####vector with cluster3 enzymes in the same order as buble plot from cluster 3
#Cluster 3 buble plot names are in SCyte_cluster3, but this table has duplicated names
#to obtain similar distribution, we remove the first entry of duplicated name
enzymesC3 <-SCyte_cluster3$Name[-which(duplicated(SCyte_cluster3$Name,fromLast = TRUE))]

#create vector with enzymes and markers order
order_enzC3 <- c(enzymesC3,markersC3)

#reorder cluster 3 rownames with order vector previously made
cluster3 <- cluster3[match(order_enzC3,rownames(cluster3)),]

annot_rowC3 <- data.frame(Enzyme=cluster3$Enzyme)
rownames(annot_rowC3) <- rownames(cluster3)
annot_rowC3$Enzyme <- factor(annot_rowC3$Enzyme, levels = c("Writer","Eraser","MK")) 

annot_colors <- list(Enzyme=c(Writer="cornflowerblue",Eraser="Red",MK="Chartreuse3"))

pheatmap(cluster3[,1:5], fontsize = 10, 
         cluster_cols = FALSE, cluster_rows = FALSE, 
         cellwidth = 30, cellheight = 15,
         annotation_row = annot_rowC3,
         annotation_colors=annot_colors,
         annotation_names_row = FALSE)

##############################################################################################
##############################################################################################
#heatmap Cluster4
cluster4 <- Enz_MK %>% filter(cell_cluster=="4")
rownames(cluster4) <-  cluster4$Name

#to sort the enzymes and MK make vector with enzymes in same order that buble plots, and add MK names

#vector with cluster 4 marker genes names
markersC4 <- rownames(cluster4[cluster4$Enzyme=="MK",])
####vector with cluster4 enzymes in the same order as buble plot from cluster 4
#Cluster 4 buble plot names are in STid_cluster4, but this table has duplicated names
#to obtain similar distribution, we remove the first entry of duplicated name
enzymesC4 <-STid_cluster4$Name[-which(duplicated(STid_cluster4$Name,fromLast = TRUE))]

#create vector with enzymes and markers order
order_enzC4 <- c(enzymesC4,markersC4)

#reorder cluster 4 rownames with order vector previously made
cluster4 <- cluster4[match(order_enzC4,rownames(cluster4)),]

annot_rowC4 <- data.frame(Enzyme=cluster4$Enzyme)
rownames(annot_rowC4) <- rownames(cluster4)
annot_rowC4$Enzyme <- factor(annot_rowC4$Enzyme, levels = c("Writer","Eraser","MK")) 

annot_colors <- list(Enzyme=c(Writer="cornflowerblue",Eraser="Red",MK="Chartreuse3"))

pheatmap(cluster4[,1:5], fontsize = 10, 
         cluster_cols = FALSE, cluster_rows = FALSE, 
         cellwidth = 30, cellheight = 15,
         annotation_row = annot_rowC4,
         annotation_colors=annot_colors,
         annotation_names_row = FALSE)

##############################################################################################
##############################################################################################
#heatmap Cluster5
cluster5 <- Enz_MK %>% filter(cell_cluster=="5")
rownames(cluster5) <-  cluster5$Name

#to sort the enzymes and MK make vector with enzymes in same order that buble plots, and add MK names

#vector with cluster 5 marker genes names
markersC5 <- rownames(cluster5[cluster5$Enzyme=="MK",])
####vector with cluster5 enzymes in the same order as buble plot from cluster 5
#Cluster 5 buble plot names are in STid_cluster5, but this table has duplicated names
#to obtain similar distribution, we remove the first entry of duplicated name
enzymesC5 <-STid_cluster5$Name[-which(duplicated(STid_cluster5$Name,fromLast = TRUE))]

#create vector with enzymes and markers order
order_enzC5 <- c(enzymesC5,markersC5)

#reorder cluster 5 rownames with order vector previously made
cluster5 <- cluster5[match(order_enzC5,rownames(cluster5)),]

annot_rowC5 <- data.frame(Enzyme=cluster5$Enzyme)
rownames(annot_rowC5) <- rownames(cluster5)
annot_rowC5$Enzyme <- factor(annot_rowC5$Enzyme, levels = c("Writer","Eraser","MK")) 

annot_colors <- list(Enzyme=c(Writer="cornflowerblue",Eraser="Red",MK="Chartreuse3"))

pheatmap(cluster5[,1:5], fontsize = 10, 
         cluster_cols = FALSE, cluster_rows = FALSE, 
         cellwidth = 30, cellheight = 15,
         annotation_row = annot_rowC5,
         annotation_colors=annot_colors,
         annotation_names_row = FALSE)
