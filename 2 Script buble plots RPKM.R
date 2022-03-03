library(tidyverse)

load("Enzymes MKs RPKM matrix.Rdata")

# Cluster tables construction for RPKM buble plots

# define histones PTMs indexes
ac_index <- c("Hac","H3ac","H4ac","H2/H4ac","H3/H4ac") # 1=H2/3/4, 3=H3, 4=H4, 24= H2/4, 34 = H3/4
me_index <- c("me1","me2","me3") #me1 me2 me3, indicates highest position described for enzyme

####################################################################################################################
####################################################################################################################

#Cluster1
cluster1 <- Enz_MK_RPKM %>% filter(cell_cluster=="1")

cluster1_Hac <- filter(cluster1,Hac %in% ac_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("Hac",sum(cluster1$Hac%in%ac_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster1_H3K4me <- filter(cluster1,H3K4me %in% me_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K4me",sum(cluster1$H3K4me%in%me_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster1_H3K36me <- filter(cluster1,H3K36me %in% me_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K36me",sum(cluster1$H3K36me%in%me_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster1_H3K9me <- filter(cluster1,H3K9me %in% me_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K9me",sum(cluster1$H3K9me%in%me_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster1_H3K27me <- filter(cluster1,H3K27me %in% me_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K27me",sum(cluster1$H3K27me%in%me_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster1_H4K20me <- filter(cluster1,H4K20me %in% me_index) %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H4K20me",sum(cluster1$H4K20me%in%me_index)), 
                                                  Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster1_H2AZK7me <- filter(cluster1,H2AZK7me=="me1") %>% 
  select(c(SGThy1,Name,Enzyme,Family)) %>% mutate(HK=rep("H2AZK7me",1), Chromatin=rep("Activation",1))

SGThy1_cluster1 <- rbind(cluster1_Hac,cluster1_H3K4me,
                         cluster1_H3K36me, cluster1_H3K9me,
                         cluster1_H3K27me, cluster1_H4K20me,
                         cluster1_H2AZK7me) 

SGThy1_cluster1$order <-  rep(nrow(SGThy1_cluster1):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
SGThy1_cluster1$HK <-  factor(SGThy1_cluster1$HK, 
                              levels = c("Hac","H3K4me","H3K36me","H3K9me","H3K27me","H4K20me","H2AZK7me"))

SGThy1_cluster1$Enzyme <- factor(SGThy1_cluster1$Enzyme, levels = c("Writer","Eraser"))

SGThy1_cluster1$Chromatin <- factor(SGThy1_cluster1$Chromatin, levels = c("Repression","Activation"))

#create plot     
ggplot(SGThy1_cluster1,
       aes(x=HK,y=reorder(Name,order),
           color=Enzyme,
           alpha=Chromatin))+
  geom_point(aes(size=SGThy1))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,5,50,100,200),
                        limits = c(.05, 300),
                        labels = c("1","5","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  scale_alpha_discrete(range=c(0.33,1))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")

####################################################################################################################
####################################################################################################################
#Cluster2
cluster2 <- Enz_MK_RPKM %>% filter(cell_cluster=="2")

cluster2_Hac <- filter(cluster2,Hac %in% ac_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")),Family) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("Hac",sum(cluster2$Hac%in%ac_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster2_H3K4me <- filter(cluster2,H3K4me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser"))) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K4me",sum(cluster2$H3K4me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster2_H3K36me <- filter(cluster2,H3K36me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K36me",sum(cluster2$H3K36me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster2_H3K9me <- filter(cluster2,H3K9me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K9me",sum(cluster2$H3K9me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster2_H3K27me <- filter(cluster2,H3K27me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K27me",sum(cluster2$H3K27me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster2_H4K20me <- filter(cluster2,H4K20me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SGKit,Name,Enzyme,Family)) %>% mutate(HK=rep("H4K20me",sum(cluster2$H4K20me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))


SGKit_cluster2 <- rbind(cluster2_Hac,cluster2_H3K4me,
                        cluster2_H3K36me, cluster2_H3K9me,
                        cluster2_H3K27me, cluster2_H4K20me)


SGKit_cluster2$order <-  rep(nrow(SGKit_cluster2):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
SGKit_cluster2$HK <-  factor(SGKit_cluster2$HK, 
                             levels = c("Hac","H3K4me","H3K36me","H3K9me","H3K27me","H4K20me"))


SGKit_cluster2$Enzyme <- factor(SGKit_cluster2$Enzyme, levels = c("Writer","Eraser"))

SGKit_cluster2$Chromatin <- factor(SGKit_cluster2$Chromatin, levels = c("Repression","Activation"))

#create plot     
ggplot(SGKit_cluster2,
       aes(x=HK,y=reorder(Name,order),
           color=Enzyme,
           alpha=Chromatin))+
  geom_point(aes(size=SGKit))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,5,50,100,200),
                        limits = c(.05, 300),
                        labels = c("1","5","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  scale_alpha_discrete(range=c(0.33,1))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")

####################################################################################################################
####################################################################################################################

#Cluster3
cluster3 <- Enz_MK_RPKM %>% filter(cell_cluster=="3")

cluster3_Hac <- filter(cluster3,Hac %in% ac_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")),Family) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("Hac",sum(cluster3$Hac%in%ac_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster3_H3K4me <- filter(cluster3,H3K4me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser"))) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K4me",sum(cluster3$H3K4me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster3_H3K36me <- filter(cluster3,H3K36me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K36me",sum(cluster3$H3K36me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster3_H3K9me <- filter(cluster3,H3K9me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K9me",sum(cluster3$H3K9me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster3_H3K27me <- filter(cluster3,H3K27me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K27me",sum(cluster3$H3K27me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster3_H4K20me <- filter(cluster3,H4K20me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(SCyte,Name,Enzyme,Family)) %>% mutate(HK=rep("H4K20me",sum(cluster3$H4K20me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))


SCyte_cluster3 <- rbind(cluster3_Hac,cluster3_H3K4me,
                        cluster3_H3K36me, cluster3_H3K9me,
                        cluster3_H3K27me, cluster3_H4K20me)


SCyte_cluster3$order <-  rep(nrow(SCyte_cluster3):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
SCyte_cluster3$HK <-  factor(SCyte_cluster3$HK, 
                             levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K9me","H3K27me","H4K20me"))

SCyte_cluster3$Enzyme <- factor(SCyte_cluster3$Enzyme, levels = c("Writer","Eraser"))

SCyte_cluster3$Chromatin <- factor(SCyte_cluster3$Chromatin, levels = c("Repression","Activation"))

#create plot     
ggplot(SCyte_cluster3,
       aes(x=HK,y=reorder(Name,order),
           color=Enzyme,
           alpha=Chromatin))+
  geom_point(aes(size=SCyte))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,5,50,100,200),
                        limits = c(.05, 300),
                        labels = c("1","5","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  scale_alpha_discrete(range=c(0.33,1))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")

####################################################################################################################
####################################################################################################################

#Cluster4
cluster4 <- Enz_MK_RPKM %>% filter(cell_cluster=="4")

cluster4_Hac <- filter(cluster4,Hac %in% ac_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")),Family) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("Hac",sum(cluster4$Hac%in%ac_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster4_H3K4me <- filter(cluster4,H3K4me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser"))) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K4me",sum(cluster4$H3K4me%in%me_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster4_H3K36me <- filter(cluster4,H3K36me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K36me",sum(cluster4$H3K36me%in%me_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster4_H3K9me <- filter(cluster4,H3K9me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K9me",sum(cluster4$H3K9me%in%me_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster4_H3K27me <- filter(cluster4,H3K27me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K27me",sum(cluster4$H3K27me%in%me_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster4_H4K20me <- filter(cluster4,H4K20me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(STid,Name,Enzyme,Family)) %>% mutate(HK=rep("H4K20me",sum(cluster4$H4K20me%in%me_index)), 
                                                Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

STid_cluster4 <- rbind(cluster4_Hac,cluster4_H3K4me,
                       cluster4_H3K36me, cluster4_H3K9me,
                       cluster4_H3K27me, cluster4_H4K20me)


STid_cluster4$order <-  rep(nrow(STid_cluster4):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
STid_cluster4$HK <-  factor(STid_cluster4$HK, 
                            levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K9me","H3K27me","H4K20me"))

STid_cluster4$Enzyme <- factor(STid_cluster4$Enzyme, levels = c("Writer","Eraser"))

STid_cluster4$Chromatin <- factor(STid_cluster4$Chromatin, levels = c("Repression","Activation"))

#create plot     
ggplot(STid_cluster4,
       aes(x=HK,y=reorder(Name,order),
           color=Enzyme,
           alpha=Chromatin))+
  geom_point(aes(size=STid))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,5,50,100,200),
                        limits = c(.05, 300),
                        labels = c("1","5","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  scale_alpha_discrete(range=c(0.33,1))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")

####################################################################################################################
####################################################################################################################

#Cluster5
cluster5 <- Enz_MK_RPKM %>% filter(cell_cluster=="5")

cluster5_Hac <- filter(cluster5,Hac %in% ac_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")),Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("Hac",sum(cluster5$Hac%in%ac_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster5_H3K4me <- filter(cluster5,H3K4me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser"))) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K4me",sum(cluster5$H3K4me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))


cluster5_H3K36me <- filter(cluster5,H3K36me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K36me",sum(cluster5$H3K36me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation", "Repression"))

cluster5_H3K9me <- filter(cluster5,H3K9me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K9me",sum(cluster5$H3K9me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster5_H3K27me <- filter(cluster5,H3K27me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K27me",sum(cluster5$H3K27me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster5_H4K20me <- filter(cluster5,H4K20me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H4K20me",sum(cluster5$H4K20me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Repression" ,"Activation"))

cluster5_H3K79me <- filter(cluster5,H3K79me %in% me_index) %>% 
  arrange(factor(Enzyme, levels = c("Writer","Eraser")), Family) %>% 
  select(c(Sperm,Name,Enzyme,Family)) %>% mutate(HK=rep("H3K79me",sum(cluster5$H4K20me%in%me_index)), 
                                                 Chromatin=ifelse(Enzyme=="Writer", "Activation" ,"Repression"))

Sperm_cluster5 <- rbind(cluster5_Hac,cluster5_H3K4me,
                        cluster5_H3K36me, cluster5_H3K9me,
                        cluster5_H3K27me, cluster5_H4K20me,
                        cluster5_H3K79me)


Sperm_cluster5$order <-  rep(nrow(Sperm_cluster5):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
Sperm_cluster5$HK <-  factor(Sperm_cluster5$HK, 
                             levels = c("Hac","H3K4me","H3K36me","H3K79me","H3K9me","H3K27me","H4K20me"))

Sperm_cluster5$Enzyme <- factor(Sperm_cluster5$Enzyme, levels = c("Writer","Eraser"))

Sperm_cluster5$Chromatin <- factor(Sperm_cluster5$Chromatin, levels = c("Repression","Activation"))

#create plot     
ggplot(Sperm_cluster5,
       aes(x=HK,y=reorder(Name,order),
           color=Enzyme,
           alpha=Chromatin))+
  geom_point(aes(size=Sperm))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,5,50,100,200),
                        limits = c(.05, 300),
                        labels = c("1","5","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  scale_alpha_discrete(range=c(0.33,1))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")