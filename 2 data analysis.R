library(DESeq2)

count.matrix <- counts[,c(1:18)]

sample.table <- as.data.frame(matrix(c("und","und","und","diff","diff","diff","prel","prel","prel",
                                         "lz","lz","lz","pd","pd","pd","rstid","rstid","rstid"), ncol = 1))
rownames(sample.table) <- colnames(count.matrix)
colnames(sample.table) <- c("group")

sample.table$batch <- rep(c("1","2","3"),6)

sample.table$group <- factor(sample.table$group, levels = c("und","diff","prel","lz","pd","rstid"))

sample.table$batch <- factor(sample.table$batch)

levels(sample.table$group)

#Check tables are OK
all(rownames(sample.table) == colnames(count.matrix))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix(countData=count.matrix, 
                                colData=sample.table, 
                                design=~batch+group)

# add analysis to dds object
dds <- DESeq(dds, test="LRT", reduced=~batch)

#see results coefficients
resultsNames(dds)

# data transformation with rlog
rld <- rlog(dds)

# see cells in PCA plot
plotPCA(rld, intgroup="group")

########################################################################################################################
# 1) SGund vs SGdiff

#obtain results
res.1 <- results(dds, contrast = c("group","und","diff"))

# convert ensembl to symbol
rownames(res.1) <- mapIds(org.Mm.eg.db, keys = rownames(res.1), keytype = "ENSEMBL", column="SYMBOL")

#inspect results fro gene of interest
res.1[which(rownames(res.1)=="Hdac2"),]

#obtain enzymes results
enz_res.1 <- res.1[which(rownames(res.1)%in%enz.genelist$Symbol),]

# upregulated in SGundiff
SGund_UP <- rownames(enz_res.1[which(enz_res.1$padj<0.05 & enz_res.1$log2FoldChange>0.5),])

# upregulated in SGdiff
SGdiff_UP <- rownames(enz_res.1[which(enz_res.1$padj<0.05 & enz_res.1$log2FoldChange< -0.5),])

#create SGund matrix for heatmap
mat.1und <- assay(rld)[which(rownames(res.1) %in% SGund_UP),1:6]
mat.1und <- mat.1und - rowMeans(mat.1und) # obtain zscore by resting the mean to each value
head(mat.1und)
rownames(mat.1und) <- mapIds(org.Mm.eg.db, keys = rownames(mat.1und), keytype = "ENSEMBL", column="SYMBOL")

# heatmap for upregulated enzymes at SGund
pheatmap(mat.1und, fontsize = 10, cellwidth = 15, cellheight = 10)

#create SGdiff matrix for heatmap
mat.1diff <- assay(rld)[which(rownames(res.1) %in% SGdiff_UP),1:6]
mat.1diff <- mat.1diff - rowMeans(mat.1diff) # obtain zscore by resting the mean to each value
rownames(mat.1diff) <- mapIds(org.Mm.eg.db, keys = rownames(mat.1diff), keytype = "ENSEMBL", column="SYMBOL")
head(mat.1diff)

# heatmap for upregulated enzymes at SGdiff
pheatmap(mat.1diff, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
# 2) Prel vs SGdiff

#obtain results
res.2 <- results(dds, contrast = c("group","prel","diff"))

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

#create PreL matrix for heatmap
mat.2UP <- assay(rld)[which(rownames(res.2) %in% Prel_UP),4:9]
mat.2UP <- mat.2UP - rowMeans(mat.2UP) # obtain zscore by resting the mean to each value
rownames(mat.2UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.2UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.2UP)

# heatmap for upregulated enzymes at PreL
pheatmap(mat.2UP, fontsize = 10, cellwidth = 15, cellheight = 10)

#dowregulated enzymes at PreL matrix
mat.2DOWN <- assay(rld)[which(rownames(res.2) %in% Prel_DOWN),4:9]
mat.2DOWN <- mat.2DOWN - rowMeans(mat.2DOWN) # obtain zscore by resting the mean to each value
rownames(mat.2DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.2DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.2DOWN)

# heatmap for downregulated enzymes at PreL
pheatmap(mat.2DOWN, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
# 3) LZ vs Prel

#obtain results
res.3 <- results(dds, contrast = c("group","lz","prel"))

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

#create upregulated LZ matrix for heatmap
mat.3UP <- assay(rld)[which(rownames(res.3) %in% LZ_UP),7:12]
mat.3UP <- mat.3UP - rowMeans(mat.3UP) # obtain zscore by resting the mean to each value
rownames(mat.3UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.3UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.3UP)

# heatmap for upregulated enzymes at LZ
pheatmap(mat.3UP, fontsize = 10,cellwidth = 15, cellheight = 10)

#create dowregulated at LZ matrix for heatmap
mat.3DOWN <- assay(rld)[which(rownames(res.3) %in% LZ_DOWN),7:12]
mat.3DOWN <- mat.3DOWN - rowMeans(mat.3DOWN) # obtain zscore by resting the mean to each value
rownames(mat.3DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.3DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.3DOWN)

# heatmap for downregulated enzymes at LZ
pheatmap(mat.3DOWN, fontsize = 10,cellwidth = 15, cellheight = 10)

########################################################################################################################
# 4) PD vs LZ

#obtain results
res.4 <- results(dds, contrast = c("group","pd","lz"))

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

#create upregulated PD matrix for heatmap
mat.4UP <- assay(rld)[which(rownames(res.4) %in% PD_UP),10:15]
mat.4UP <- mat.4UP - rowMeans(mat.4UP) # obtain zscore by resting the mean to each value
rownames(mat.4UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.4UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.4UP)

# heatmap for upregulated enzymes at PD
pheatmap(mat.4UP, fontsize = 10, cellwidth = 15, cellheight = 10)

#create downregulated PD matrix for heatmap
mat.4DOWN <- assay(rld)[which(rownames(res.4) %in% PD_DOWN),10:15]
mat.4DOWN <- mat.4DOWN - rowMeans(mat.4DOWN) # obtain zscore by resting the mean to each value
rownames(mat.4DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.4DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.4DOWN)

# heatmap for downregulated enzymes at PD
pheatmap(mat.4DOWN, fontsize = 10, cellwidth = 15, cellheight = 10)

########################################################################################################################
# 5) PD vs RStid

res.5 <- results(dds, contrast = c("group","rstid","pd"))

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

#create upregulated RStid matrix for heatmap
mat.5UP <- assay(rld)[which(rownames(res.5) %in% RStid_UP),13:18]
mat.5UP <- mat.5UP - rowMeans(mat.5UP) # obtain zscore by resting the mean to each value
rownames(mat.5UP) <- mapIds(org.Mm.eg.db, keys = rownames(mat.5UP), keytype = "ENSEMBL", column="SYMBOL")
head(mat.5UP)

# heatmap for upregulated enzymes at RStid
pheatmap(mat.5UP, fontsize = 10,cellwidth = 15, cellheight = 10)

#create downregulated RStid matrix for heatmap
mat.5DOWN <- assay(rld)[which(rownames(res.5) %in% RStid_DOWN),13:18]
mat.5DOWN <- mat.5DOWN - rowMeans(mat.5DOWN) # obtain zscore by resting the mean to each value
rownames(mat.5DOWN) <- mapIds(org.Mm.eg.db, keys = rownames(mat.5DOWN), keytype = "ENSEMBL", column="SYMBOL")
head(mat.5DOWN)

# heatmap for downregulated enzymes at RStid
pheatmap(mat.5DOWN, fontsize = 10,cellwidth = 15, cellheight = 10)

########################################################################################################################
"%notin%" <- Negate("%in%")

#search for enzymes not upregulated at any stage
unassigned_enz <- enz.genelist$Symbol[which(enz.genelist$Symbol%notin%SGund_UP &
                                          enz.genelist$Symbol%notin%SGdiff_UP&
                                          enz.genelist$Symbol%notin%Prel_UP&
                                          enz.genelist$Symbol%notin%LZ_UP&
                                          enz.genelist$Symbol%notin%PD_UP&
                                          enz.genelist$Symbol%notin%RStid_UP)]

unassigned_enz_counts <- counts[which(counts$symbol%in%unassigned_enz),]
rownames(unassigned_enz_counts) <- unassigned_enz_counts$symbol

pheatmap(unassigned_enz_counts[,1:18], fontsize = 10, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE)
