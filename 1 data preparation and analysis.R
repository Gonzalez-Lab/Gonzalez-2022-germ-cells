library(tidyverse)
library(org.Mm.eg.db)
library(pheatmap)
library(DESeq2)

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
gene <- "Kmt2d"
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

########################################################################################################################
# 1) SGund vs SGdiff
count.matrix.1 <- counts[,c(1,2,3,4,5,6)]

sample.table.1 <- as.data.frame(matrix(c("und","und","und","diff","diff","diff"), ncol = 1))
rownames(sample.table.1) <- colnames(count.matrix.1)
colnames(sample.table.1) <- c("group")

sample.table.1$group <- factor(sample.table.1$group)

levels(sample.table.1$group)

#Check tables are OK
all(rownames(sample.table.1) == colnames(count.matrix.1))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds.1 <- DESeqDataSetFromMatrix(countData=count.matrix.1, 
                              colData=sample.table.1, 
                              design=~group)

# add analysis to dds.1 object
dds.1 <- estimateSizeFactors(dds.1)
dds.1 <- DESeq(dds.1)

#obtain results
res.1 <- results(dds.1)

table(res.1$padj != "NA")
table(res.1$padj < 0.05)

# convert ensembl to symbol
rownames(res.1) <- mapIds(org.Mm.eg.db, keys = rownames(res.1), keytype = "ENSEMBL", column="SYMBOL")

#inspect results fro gene of interes
res.1[which(rownames(res.1)=="Mcm3ap"),]

#obtain enzymes results
enz_res.1 <- res.1[which(rownames(res.1)%in%enz.genelist$Symbol),]

# upregulated in SGundiff
SGund_UP <- rownames(enz_res.1[which(enz_res.1$padj<0.05 & enz_res.1$log2FoldChange>0.5),])

# upregulated in SGdiff
SGdiff_UP <- rownames(enz_res.1[which(enz_res.1$padj<0.05 & enz_res.1$log2FoldChange< -0.5),])

# data transformation with rlog
rld.1 <- rlog(dds.1)

plotPCA(rld.1, intgroup="group")

#create SGund matrix for heatmap
mat.1und <- assay(rld.1)[which(rownames(res.1) %in% c(SGund_UP)),]
mat.1und <- mat.1und - rowMeans(mat.1und) # obtain zscore by resting the mean to each value
head(mat.1und)
rownames(mat.1und) <- mapIds(org.Mm.eg.db, keys = rownames(mat.1und), keytype = "ENSEMBL", column="SYMBOL")

# heatmap for upregulated enzymes at SGund
pheatmap(mat.1und, fontsize = 10, cellwidth = 15, cellheight = 10)

#create SGdiff matrix for heatmap
mat.1diff <- assay(rld.1)[which(rownames(res.1) %in% c(SGdiff_UP)),]
mat.1diff <- mat.1diff - rowMeans(mat.1diff) # obtain zscore by resting the mean to each value
rownames(mat.1diff) <- mapIds(org.Mm.eg.db, keys = rownames(mat.1diff), keytype = "ENSEMBL", column="SYMBOL")
head(mat.1diff)

# heatmap for upregulated enzymes at SGdiff
pheatmap(mat.1diff, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
########################################################################################################################
# 2) Prel vs SGdiff
count.matrix.2 <- counts[,c(4,5,6,7,8,9)]

sample.table.2 <- as.data.frame(matrix(c("diff","diff","diff","prel","prel","prel"), ncol = 1))
rownames(sample.table.2) <- colnames(count.matrix.2)
colnames(sample.table.2) <- c("group")

sample.table.2$group <- factor(sample.table.2$group)

levels(sample.table.2$group)

#check tables are OK
all(rownames(sample.table.2) == colnames(count.matrix.2))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds.2 <- DESeqDataSetFromMatrix(countData=count.matrix.2, 
                                colData=sample.table.2, 
                                design=~group)

# add analysis to dds.2 object
dds.2 <- estimateSizeFactors(dds.2)
dds.2 <- DESeq(dds.2)

#obtain results
res.2 <- results(dds.2)

table(res.2$padj != "NA")
table(res.2$padj < 0.05)

#change ensembl to symbol
rownames(res.2) <- mapIds(org.Mm.eg.db, keys = rownames(res.2), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interes
res.2[which(rownames(res.2)=="Kat6b"),]

#obtain enzymes results
enz_res.2 <- res.2[which(rownames(res.2)%in%enz.genelist$Symbol),]

# upregulated in Prel
Prel_UP <- rownames(enz_res.2[which(enz_res.2$padj<0.05 & enz_res.2$log2FoldChange>0.5),])

# upregulated in SGdiff / downregulated in PreL
Prel_DOWN <- rownames(enz_res.2[which(enz_res.2$padj<0.05 & enz_res.2$log2FoldChange< -0.5),])

# data transformation with rlog
rld.2 <- rlog(dds.2)

plotPCA(rld.2, intgroup="group")

#create PreL matrix for heatmap
mat.2UP <- assay(rld.2)[which(rownames(res.2) %in% c(Prel_UP)),]
mat.2UP <- mat.2UP - rowMeans(mat.2UP) # obtain zscore by resting the mean to each value
rownames(mat.2UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.2UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.2UP)

# heatmap for upregulated enzymes at PreL
pheatmap(mat.2UP, fontsize = 10, cellwidth = 15, cellheight = 10)

#dowregulated enzymes at PreL matrix
mat.2DOWN <- assay(rld.2)[which(rownames(res.2) %in% c(Prel_DOWN)),]
mat.2DOWN <- mat.2DOWN - rowMeans(mat.2DOWN) # obtain zscore by resting the mean to each value
rownames(mat.2DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.2DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.2DOWN)

# heatmap for downregulated enzymes at PreL
pheatmap(mat.2DOWN, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
########################################################################################################################
# 3) LZ vs Prel
count.matrix.3 <- counts[,c(7,8,9,10,11,12)]

sample.table.3 <- as.data.frame(matrix(c("prel","prel","prel","LZ","LZ","LZ"), ncol = 1))
rownames(sample.table.3) <- colnames(count.matrix.3)
colnames(sample.table.3) <- c("group")

sample.table.3$group <- factor(sample.table.3$group)

levels(sample.table.3$group)
sample.table.3$group <- relevel(sample.table.3$group, "prel")

#check tables are OK
all(rownames(sample.table.3) == colnames(count.matrix.3))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds.3 <- DESeqDataSetFromMatrix(countData=count.matrix.3, 
                                colData=sample.table.3, 
                                design=~group)

# add analysis to dds.3 object
dds.3 <- estimateSizeFactors(dds.3)
dds.3 <- DESeq(dds.3)

#obtain results
res.3 <- results(dds.3)

table(res.3$padj != "NA")
table(res.3$padj < 0.05)

#change ensembl to symbol
rownames(res.3) <- mapIds(org.Mm.eg.db, keys = rownames(res.3), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
res.3[which(rownames(res.3)=="Nsd2"),]

#obtain enzymes results
enz_res.3 <- res.3[which(rownames(res.3)%in%enz.genelist$Symbol),]

# upregulated in LZ
LZ_UP <- rownames(enz_res.3[which(enz_res.3$padj<0.05 & enz_res.3$log2FoldChange>0.5),])

# upregulated in PreL / Downregulated in LZ
LZ_DOWN <- rownames(enz_res.3[which(enz_res.3$padj<0.05 & enz_res.3$log2FoldChange< -0.5),])

# data transformation with rlog
rld.3 <- rlog(dds.3)

plotPCA(rld.3, intgroup="group")

#create upregulated LZ matrix for heatmap
mat.3UP <- assay(rld.3)[which(rownames(res.3) %in% c(LZ_UP)),]
mat.3UP <- mat.3UP - rowMeans(mat.3UP) # obtain zscore by resting the mean to each value
rownames(mat.3UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.3UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.3UP)

# heatmap for upregulated enzymes at LZ
pheatmap(mat.3UP, fontsize = 10,cellwidth = 15, cellheight = 10)

#create dowregulated at LZ matrix for heatmap
mat.3DOWN <- assay(rld.3)[which(rownames(res.3) %in% c(LZ_DOWN)),]
mat.3DOWN <- mat.3DOWN - rowMeans(mat.3DOWN) # obtain zscore by resting the mean to each value
rownames(mat.3DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.3DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.3DOWN)

# heatmap for downregulated enzymes at LZ
pheatmap(mat.3DOWN, fontsize = 10,cellwidth = 15, cellheight = 10)

########################################################################################################################
########################################################################################################################
# 4) PD vs LZ
count.matrix.4 <- counts[,c(10,11,12,13,14,15)]

sample.table.4 <- as.data.frame(matrix(c("LZ","LZ","LZ","PD","PD","PD"), ncol = 1))
rownames(sample.table.4) <- colnames(count.matrix.4)
colnames(sample.table.4) <- c("group")

sample.table.4$group <- factor(sample.table.4$group)

levels(sample.table.4$group)

#check tables are OK
all(rownames(sample.table.4) == colnames(count.matrix.4))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds.4 <- DESeqDataSetFromMatrix(countData=count.matrix.4, 
                                colData=sample.table.4, 
                                design=~group)

# add analysis to dds.4 object
dds.4 <- estimateSizeFactors(dds.4)
dds.4 <- DESeq(dds.4)

#obtain results
res.4 <- results(dds.4)

table(res.4$padj != "NA")
table(res.4$padj < 0.05)

#change ensembl to symbol
rownames(res.4) <- mapIds(org.Mm.eg.db, keys = rownames(res.4), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
res.4[which(rownames(res.4)=="Hdac3"),]

#obtain enzymes results
enz_res.4 <- res.4[which(rownames(res.4)%in%enz.genelist$Symbol),]

# upregulated in PD
PD_UP <- rownames(enz_res.4[which(enz_res.4$padj<0.05 & enz_res.4$log2FoldChange>0.5),])

# upregulated in LZ / Downregulated in PD
PD_DOWN <- rownames(enz_res.4[which(enz_res.4$padj<0.05 & enz_res.4$log2FoldChange< -0.5),])

# data transformation with rlog
rld.4 <- rlog(dds.4)

plotPCA(rld.4, intgroup="group")

#create upregulated PD matrix for heatmap
mat.4UP <- assay(rld.4)[which(rownames(res.4) %in% c(PD_UP)),]
mat.4UP <- mat.4UP - rowMeans(mat.4UP) # obtain zscore by resting the mean to each value
rownames(mat.4UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.4UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.4UP)

# heatmap for upregulated enzymes at PD
pheatmap(mat.4UP, fontsize = 10, cellwidth = 15, cellheight = 10)

#create downregulated PD matrix for heatmap
mat.4DOWN <- assay(rld.4)[which(rownames(res.4) %in% c(PD_DOWN)),]
mat.4DOWN <- mat.4DOWN - rowMeans(mat.4DOWN) # obtain zscore by resting the mean to each value
rownames(mat.4DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.4DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.4DOWN)

# heatmap for downregulated enzymes at PD
pheatmap(mat.4DOWN, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
########################################################################################################################
# 5) PD vs RStid
count.matrix.5 <- counts[,c(13,14,15,16,17,18)]

sample.table.5 <- as.data.frame(matrix(c("PD","PD","PD","rstid","rstid","rstid"), ncol = 1))
rownames(sample.table.5) <- colnames(count.matrix.5)
colnames(sample.table.5) <- c("group")

sample.table.5$group <- factor(sample.table.5$group)

levels(sample.table.5$group)

#check tables are OK
all(rownames(sample.table.5) == colnames(count.matrix.5))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds.5 <- DESeqDataSetFromMatrix(countData=count.matrix.5, 
                                colData=sample.table.5, 
                                design=~group)

# add analysis to dds.5 object
dds.5 <- estimateSizeFactors(dds.5)
dds.5 <- DESeq(dds.5)

#obtain results
res.5 <- results(dds.5)

table(res.5$padj != "NA")
table(res.5$padj < 0.05)

#change ensembl to symbol
rownames(res.5) <- mapIds(org.Mm.eg.db, keys = rownames(res.5), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
res.5[which(rownames(res.5)=="Hdac3"),]

#obtain enzymes results
enz_res.5 <- res.5[which(rownames(res.5)%in%enz.genelist$Symbol),]

# upregulated in RStid
RStid_UP <- rownames(enz_res.5[which(enz_res.5$padj<0.05 & enz_res.5$log2FoldChange>0.5),])

# upregulated in PD / Downregulated in RStid
RStid_DOWN <- rownames(enz_res.5[which(enz_res.5$padj<0.05 & enz_res.5$log2FoldChange< -0.5),])

# data transformation with rlog
rld.5 <- rlog(dds.5)

plotPCA(rld.5, intgroup="group")

#create upregulated RStid matrix for heatmap
mat.5UP <- assay(rld.5)[which(rownames(res.5) %in% c(RStid_UP)),]
mat.5UP <- mat.5UP - rowMeans(mat.5UP) # obtain zscore by resting the mean to each value
rownames(mat.5UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.5UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.5UP)

# heatmap for upregulated enzymes at RStid
pheatmap(mat.5UP, fontsize = 10,cellwidth = 15, cellheight = 10)

#create downregulated RStid matrix for heatmap
mat.5DOWN <- assay(rld.5)[which(rownames(res.5) %in% c(RStid_DOWN)),]
mat.5DOWN <- mat.5DOWN - rowMeans(mat.5DOWN) # obtain zscore by resting the mean to each value
rownames(mat.5DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.5DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.5DOWN)

# heatmap for downregulated enzymes at RStid
pheatmap(mat.5DOWN, fontsize = 10,cellwidth = 15, cellheight = 10)

########################################################################################################################
"%notin%" <- Negate("%in%")

#search for enzymes not upregulated at any stage
unassigned_enz <- genelist$Symbol[which(genelist$Symbol%notin%SGund_UP &
                                          genelist$Symbol%notin%SGdiff_UP&
                                          genelist$Symbol%notin%Prel_UP&
                                          genelist$Symbol%notin%LZ_UP&
                                          genelist$Symbol%notin%PD_UP&
                                          genelist$Symbol%notin%RStid_UP)]

unassigned_enz_counts <- counts[which(counts$symbol%in%unassigned_enz),]
rownames(unassigned_enz_counts) <- unassigned_enz_counts$symbol

pheatmap(unassigned_enz_counts[,1:18], fontsize = 10, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE)
