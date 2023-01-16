library(tidyverse)
library(org.Mm.eg.db)
library(pheatmap)

#upload the count matrix
counts <-  read.csv("GSE162740 germ cells counts.csv", header=TRUE)

counts$symbol <- mapIds(org.Mm.eg.db, keys = counts$Gene, keytype = "ENSEMBL", column="SYMBOL")

#upload gene lengths
gene.lengths <- read.csv("gencodevM29.genelength.csv", header=TRUE)

counts$gene.lengths <- gene.lengths$mean[match(counts$Gene,gene.lengths$gene)]

counts <- counts %>% na.omit()
rownames(counts)  <- counts$Gene
counts <- counts[,-1]

#look for gene of interest
gene <- "Brdt"
counts[counts$symbol==gene,]

# upload list for epigenetic enzymes genes
enz.genelist <- read.csv("epigenetic enzymes genes mm10.csv", header=TRUE)

# upload list for MK genes
MK.genelist <- read.csv("MK genes mm10.csv", header=TRUE)

# upload list for somatic cells MK genes
somaticMK.genelist <- read.csv("somatic MK genelist mm10.csv", header=TRUE)

# counts to RPKM conversion
RPKMs <- as.data.frame((counts[,1]/(sum(counts[,1])/1000000))/(counts$gene.lengths/1000))
colnames(RPKMs) <- "SGund.1"
rownames(RPKMs) <- rownames(counts)
RPKMs$SGund.2 <- (counts[,2]/(sum(counts[,2])/1000000))/(counts$gene.lengths/1000)
RPKMs$SGund.3 <- (counts[,3]/(sum(counts[,3])/1000000))/(counts$gene.lengths/1000)
RPKMs$SGdiff.1 <- (counts[,4]/(sum(counts[,4])/1000000))/(counts$gene.lengths/1000)
RPKMs$SGdiff.2 <- (counts[,5]/(sum(counts[,5])/1000000))/(counts$gene.lengths/1000)
RPKMs$SGdiff.3 <- (counts[,6]/(sum(counts[,6])/1000000))/(counts$gene.lengths/1000)
RPKMs$PreL.1 <- (counts[,7]/(sum(counts[,7])/1000000))/(counts$gene.lengths/1000)
RPKMs$PreL.2 <- (counts[,8]/(sum(counts[,8])/1000000))/(counts$gene.lengths/1000)
RPKMs$PreL.3 <- (counts[,9]/(sum(counts[,9])/1000000))/(counts$gene.lengths/1000)
RPKMs$LZ.1 <- (counts[,10]/(sum(counts[,10])/1000000))/(counts$gene.lengths/1000)
RPKMs$LZ.2 <- (counts[,11]/(sum(counts[,11])/1000000))/(counts$gene.lengths/1000)
RPKMs$LZ.3 <- (counts[,12]/(sum(counts[,12])/1000000))/(counts$gene.lengths/1000)
RPKMs$PD.1 <- (counts[,13]/(sum(counts[,13])/1000000))/(counts$gene.lengths/1000)
RPKMs$PD.2 <- (counts[,14]/(sum(counts[,14])/1000000))/(counts$gene.lengths/1000)
RPKMs$PD.3 <- (counts[,15]/(sum(counts[,15])/1000000))/(counts$gene.lengths/1000)
RPKMs$RStid.1 <- (counts[,16]/(sum(counts[,16])/1000000))/(counts$gene.lengths/1000)
RPKMs$RStid.2 <- (counts[,17]/(sum(counts[,17])/1000000))/(counts$gene.lengths/1000)
RPKMs$RStid.3 <- (counts[,18]/(sum(counts[,18])/1000000))/(counts$gene.lengths/1000)
RPKMs$symbol <- counts$symbol

RPKMs <- RPKMs %>% distinct(symbol, .keep_all= TRUE)
rownames(RPKMs) <- RPKMs$symbol
RPKMs <- RPKMs[,-19]

# epigenetic enzymes genes RPKMs table
enz_RPKM <- RPKMs[which(rownames(RPKMs) %in% enz.genelist$Symbol),]

# heatmap for epigenetic enzymes zscore
pheatmap(enz_RPKM, scale = "row", fontsize_row = 7)

# Developmental marker genes RPKMs table
MK_RPKM <- RPKMs[which(rownames(RPKMs) %in% MK.genelist$Symbol),]
MK_RPKM <- MK_RPKM[match(MK.genelist$Symbol,rownames(MK_RPKM)),]

# heatmap for developmental markers zscore
pheatmap(MK_RPKM, scale = "row", fontsize_row = 7, cutree_rows = 1)

# Marker genes for somatic cells
Leydig_Sertoli_RPKM <- RPKMs[which(rownames(RPKMs) %in% somaticMK.genelist$symbol),]

