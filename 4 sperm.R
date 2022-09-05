library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggwordcloud)
library(org.Mm.eg.db)

#upload tables
enz.genelist <- read.csv("epigenetic enzymes genes mm10.csv", header=TRUE)
gene.lengths <- read.csv("gencodevM29.genelength.csv", header=TRUE)
gene.lengths$symbol <- mapIds(org.Mm.eg.db, keys = gene.lengths$gene, keytype = "ENSEMBL", column="SYMBOL")

##########################################################################################################
#table 1 is total and head sperm mRNA from GSE81216
sperm_table1 <- read.delim("GSE81216 sperm counts.tab", header=TRUE)

sperm_table1$gene.lengths <- gene.lengths$mean[match(sperm_table1$GeneId,gene.lengths$symbol)]

table1_rpkm <- cbind(total.1=(sperm_table1[,2]/(sum(sperm_table1[,2])/1000000))/(sperm_table1$gene.lengths/1000),
                     total.2=(sperm_table1[,3]/(sum(sperm_table1[,3])/1000000))/(sperm_table1$gene.lengths/1000),
                     head.1=(sperm_table1[,4]/(sum(sperm_table1[,4])/1000000))/(sperm_table1$gene.lengths/1000),
                     head.2=(sperm_table1[,5]/(sum(sperm_table1[,5])/1000000))/(sperm_table1$gene.lengths/1000))

table1_rpkm <- as.data.frame(table1_rpkm) 
table1_rpkm$mean_rpkm <- rowMeans(table1_rpkm)
table1_rpkm$symbol <- sperm_table1$GeneId

enz_table1_rpkm <- table1_rpkm[which(table1_rpkm$symbol %in% enz.genelist$Symbol),]

enz_table1_rpkm <- enz_table1_rpkm[match(enz_category$symbol, enz_table1_rpkm$symbol),] 

enz_table1_rpkm$category <- enz_category$pos_trasl

enz_table1_rpkm$percent_rank <- percent_rank(enz_table1_rpkm$mean_rpkm)

###########################################################################################################
##########################################################################################################
#table 2 is total sperm mRNA from GSE88732
sperm_table2 <- read.delim("GSE88732 sperm counts.tab", header=TRUE)

sperm_table2$gene.lengths <- gene.lengths$mean[match(sperm_table2$GeneId,gene.lengths$symbol)]

table2_rpkm <- cbind(total.1=(sperm_table2[,2]/(sum(sperm_table2[,2])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.2=(sperm_table2[,3]/(sum(sperm_table2[,3])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.3=(sperm_table2[,4]/(sum(sperm_table2[,4])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.4=(sperm_table2[,5]/(sum(sperm_table2[,5])/1000000))/(sperm_table2$gene.lengths/1000))


table2_rpkm <- as.data.frame(table2_rpkm) 
table2_rpkm$mean_rpkm <- rowMeans(table2_rpkm)
table2_rpkm$symbol <- sperm_table2$GeneId

enz_table2_rpkm <- table2_rpkm[which(table2_rpkm$symbol %in% enz.genelist$Symbol),]

enz_table2_rpkm <- enz_table2_rpkm[match(enz_category$symbol, enz_table2_rpkm$symbol),] 

enz_table2_rpkm$category <- enz_category$pos_trasl

enz_table2_rpkm$percent_rank <- percent_rank(enz_table2_rpkm$mean_rpkm)

###########################################################################################################
##########################################################################################################
#table 3 is total sperm mRNA from E-MTAB-5834 control vs MSUS mice
sperm_table3 <- read.delim("E-MTAB-5834 sperm counts.tab", header=TRUE)

sperm_table3$gene.lengths <- gene.lengths$mean[match(sperm_table3$GeneId,gene.lengths$symbol)]

table3_rpkm <- cbind(total.1=(sperm_table3[,2]/(sum(sperm_table3[,2])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.2=(sperm_table3[,3]/(sum(sperm_table3[,3])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.3=(sperm_table3[,4]/(sum(sperm_table3[,4])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.4=(sperm_table3[,5]/(sum(sperm_table3[,5])/1000000))/(sperm_table3$gene.lengths/1000))


table3_rpkm <- as.data.frame(table3_rpkm) 
table3_rpkm$mean_rpkm <- rowMeans(table3_rpkm)
table3_rpkm$symbol <- sperm_table3$GeneId

enz_table3_rpkm <- table3_rpkm[which(table3_rpkm$symbol %in% enz.genelist$Symbol),]

enz_table3_rpkm <- enz_table3_rpkm[match(enz_category$symbol, enz_table3_rpkm$symbol),] 

enz_table3_rpkm$category <- enz_category$pos_trasl

enz_table3_rpkm$percent_rank <- percent_rank(enz_table3_rpkm$mean_rpkm)

###############################################################################################################
#mean percent rank for the three studies
sperm_tables <- as.data.frame(cbind(table1=enz_table1_rpkm$percent_rank, 
                                    table2=enz_table2_rpkm$percent_rank,
                                    table3=enz_table3_rpkm$percent_rank))

sperm_tables$mean_percent_rank <- rowMeans(sperm_tables)

sperm_tables$symbol <- enz_table1_rpkm$symbol

sperm_tables$category <- enz_table1_rpkm$category

sperm_tables <- sperm_tables[order(sperm_tables$mean_percent_rank),]

ggplot(data = sperm_tables, 
       aes(label = symbol, size = mean_percent_rank, col=as.character(category), alpha=mean_percent_rank)) + 
  geom_text_wordcloud(rm_outside = TRUE, max_steps = 2,
                      grid_size = 5, eccentricity =0.9)+
  theme_void()+
  scale_color_manual(values = c("black", "#d52954"))

###############################################################################################################
# correlation sperm vs RStid

stid <- percent_rank(rowMeans(enz_RPKM[,16:18]))
stid <- stid[match(sperm_tables$symbol,names(stid))]

stid_sperm <- as.data.frame(cbind(stid,sperm=sperm_tables$mean_percent_rank))

ggscatter(stid_sperm, 
          x = "stid", 
          y = "sperm", 
          add = "reg.line", 
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "mean percent rank RStid", 
          ylab = "mean percent rank sperm",
          title = "",
          color = "black") + theme_bw() + theme(text = element_text(size = 10))


