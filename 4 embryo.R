library(tidyverse)
library(org.Mm.eg.db)

#Ribosome-bound RNA at 1 cell embryo
cell1rep1 <- read.delim("liRiboseq_1-cell_rep1.quant.genes.txt", header=TRUE)
cell1rep2 <- read.delim("liRiboseq_1-cell_rep2.quant.genes.txt", header=TRUE)

#create 1 cell embryo riboseq table
embryo_1c <- cbind(rep1=cell1rep1$FPKM, rep2=cell1rep2$FPKM)
embryo_1c <- cbind(embryo_1c, mean=rowMeans(embryo_1c))
embryo_1c <- as.data.frame(cbind(embryo_1c,ensembl=cell1rep1$gene_id))

embryo_1c$symbol <- mapIds(org.Mm.eg.db, keys = substr(embryo_1c$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
embryo_1c[which(embryo_1c$symbol=="N6amt1"),]

#total RNA at 1 cell embryo
cell1rep1T <- read.delim("totalRNA_1-cell_rep1.quant.genes.txt", header=TRUE)
cell1rep2T <- read.delim("totalRNA_1-cell_rep2.quant.genes.txt", header=TRUE)

#create 1 cell embryo total RNA table
embryo_1cT <- cbind(rep1=cell1rep1T$FPKM, rep2=cell1rep2T$FPKM)
embryo_1cT <- cbind(embryo_1cT, mean=rowMeans(embryo_1cT))
embryo_1cT <- as.data.frame(cbind(embryo_1cT,ensembl=cell1rep1T$gene_id))

embryo_1cT$symbol <- mapIds(org.Mm.eg.db, keys = substr(embryo_1cT$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
embryo_1cT[which(embryo_1cT$symbol=="N6amt1"),]

embryo_1cT <- embryo_1cT[which(embryo_1cT$ensembl %in% embryo_1c$ensembl),]

embryo_1cT <- embryo_1cT[match(embryo_1c$symbol, embryo_1cT$symbol),]

RPKM_TE_table <- cbind(totalRNA=embryo_1cT$mean,translRNA=embryo_1c$mean)
RPKM_TE_table <- cbind(RPKM_TE_table, TE=as.numeric(RPKM_TE_table[,2])/as.numeric(RPKM_TE_table[,1]))
RPKM_TE_table <- as.data.frame(RPKM_TE_table)
RPKM_TE_table$symbol <- embryo_1cT$symbol

RPKM_TE_table[which(RPKM_TE_table$symbol=="N6amt1"),]

#obtain enzymes results
enz_RPKM_TE_table <- RPKM_TE_table[which(RPKM_TE_table$symbol %in% enz.genelist$Symbol),]

#create enzymes category table for positive translation (TE>0) vs not translated (TE=0) 
enz_category <- as.data.frame(cbind(symbol=enz_RPKM_TE_table$symbol, pos_trasl=ifelse(enz_RPKM_TE_table$TE>0,1,0)))


########################################################################################################
####################################################################################################
#MII embryo analysis

#total RNA at MII oocyte
MIIrep1T <- read.delim("totalRNA_mII_rep1.quant.genes.txt", header=TRUE)
MIIrep2T <- read.delim("totalRNA_MII_rep2.quant.genes.txt", header=TRUE)
#
MII <- as.data.frame(cbind(mii.rep1=MIIrep1T$FPKM, mii.rep2=MIIrep2T$FPKM, ensembl=MIIrep1T$gene_id))
MII$symbol <- mapIds(org.Mm.eg.db, keys = substr(MII$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

enz.MII <- MII[which(MII$symbol %in% enz.genelist$Symbol),]

enz.embryo_1cT <- embryo_1cT[which(embryo_1cT$symbol %in% enz.genelist$Symbol),]

enz.MII <- enz.MII[which(enz.MII$symbol%in%enz.embryo_1cT$symbol),]
enz.MII <- enz.MII[match(enz.MII$symbol,enz.embryo_1cT$symbol),]

MII.embryo <- cbind(mii.rep1 = as.numeric(enz.MII$mii.rep1),
                    Mii.rep2 = as.numeric(enz.MII$mii.rep2),
                    embryo.rep1 = as.numeric(enz.embryo_1cT$rep1),
                    embryo.rep2 = as.numeric(enz.embryo_1cT$rep2))

MII.embryo <- cbind(MII.embryo,
                    meanMII = rowMeans(MII.embryo[,1:2]))

MII.embryo <- as.data.frame(MII.embryo)
MII.embryo$symbol <- enz.MII$symbol

MII.embryo.stat <- cbind(embryo_1 = MII.embryo[,3]/MII.embryo[,5],
                         embryo_2 = MII.embryo[,4]/MII.embryo[,5])

MII.embryo.stat <- as.data.frame(MII.embryo.stat)
rownames(MII.embryo.stat) <- enz.MII$symbol 
MII.embryo.stat <- MII.embryo.stat[-which(rownames(MII.embryo.stat)=="Prdm8"),]

enzyme <- "Rsbn1"

t.test(MII.embryo.stat[which(rownames(MII.embryo.stat)==enzyme),],mu=1, alternative = "greater")$p.value
t.test(MII.embryo.stat[which(rownames(MII.embryo.stat)==enzyme),],mu=1, alternative = "greater")$stat