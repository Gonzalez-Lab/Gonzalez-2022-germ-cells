library(tidyverse)
library(ggrepel)

# define histones PTMs indexes
ac_index <- c("Hac","H3ac","H4ac","H2/H4ac","H3/H4ac") 
me_index <- c("me1","me2","me3") #me1 me2 me3, indicates highest position described for enzyme

####################################################################################################################
####################################################################################################################
#the enz.genelist object contains the enzymes metadata

#SGund enzymes: filter SGund enzymes from enz.genelist file and add the mean RPKM value at SGund
Enz_SGund <- enz.genelist[which(enz.genelist$Symbol %in% SGund_UP),]
Enz_SGund_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% SGund_UP,1:3])
Enz_SGund$RPKM <- Enz_SGund_RPKM[match(Enz_SGund$Symbol,names(Enz_SGund_RPKM))] 

#construct table for dotplot for each mark
Enz_SGund_Hac <- filter(Enz_SGund,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme, Hac)) %>% mutate(HK=rep("Hac",sum(Enz_SGund$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_SGund_H3K4me <- filter(Enz_SGund,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme, H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_SGund$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_SGund_H3K36me <- filter(Enz_SGund,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme, H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_SGund$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_SGund_H3K9me <- filter(Enz_SGund,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_SGund$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_SGund_H3K27me <- filter(Enz_SGund,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_SGund$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_SGund_H4K20me <- filter(Enz_SGund,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_SGund$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_SGund_H3K79me <- filter(Enz_SGund,H3K79me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_SGund$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_SGund_H2AZK7me <- filter(Enz_SGund,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_SGund_H4K12me <- filter(Enz_SGund,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")


SGund <- rbind(Enz_SGund_Hac,Enz_SGund_H3K4me,
              Enz_SGund_H3K36me, Enz_SGund_H3K9me,
              Enz_SGund_H3K27me, Enz_SGund_H4K20me,
              Enz_SGund_H3K79me, Enz_SGund_H2AZK7me,
              Enz_SGund_H4K12me) 

#create column with descending order to arrange graph
SGund$order <-  rep(nrow(SGund):1) 

#Factorize variables so they appear in the desired order
SGund$HK <-  factor(SGund$HK, 
                    levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

SGund$Enzyme <- factor(SGund$Enzyme, levels = c("Writer","Eraser"))

#create plot SGund  
ggplot(SGund,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
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
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)


####################################################################################################################
####################################################################################################################

#SGdiff enzymes: filter SGdiff enzymes from enz.genelist file and add the mean RPKM value at SGdiff
Enz_SGdiff <- enz.genelist[which(enz.genelist$Symbol %in% SGdiff_UP),]
Enz_SGdiff_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% SGdiff_UP,4:6])
Enz_SGdiff$RPKM <- Enz_SGdiff_RPKM[match(Enz_SGdiff$Symbol,names(Enz_SGdiff_RPKM))] 

#construct table for dotplot for each mark
Enz_SGdiff_Hac <- filter(Enz_SGdiff,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_SGdiff$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_SGdiff_H3K4me <- filter(Enz_SGdiff,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme, H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_SGdiff$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_SGdiff_H3K36me <- filter(Enz_SGdiff,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_SGdiff$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_SGdiff_H3K9me <- filter(Enz_SGdiff,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_SGdiff$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_SGdiff_H3K27me <- filter(Enz_SGdiff,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_SGdiff$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_SGdiff_H4K20me <- filter(Enz_SGdiff,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_SGdiff$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_SGdiff_H3K79me <- filter(Enz_SGdiff,H3K79me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Family,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_SGdiff$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_SGdiff_H2AZK7me <- filter(Enz_SGdiff,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_SGdiff_H4K12me <- filter(Enz_SGdiff,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

SGdiff <- rbind(Enz_SGdiff_Hac,Enz_SGdiff_H3K4me,
               Enz_SGdiff_H3K36me, Enz_SGdiff_H3K9me,
               Enz_SGdiff_H3K27me, Enz_SGdiff_H4K20me,
               Enz_SGdiff_H3K79me, Enz_SGdiff_H2AZK7me,
               Enz_SGdiff_H4K12me) 

#create column with descending order to arrange graph
SGdiff$order <-  rep(nrow(SGdiff):1) 

#Factorize variables so they appear in the desired order
SGdiff$HK <-  factor(SGdiff$HK, 
                    levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

SGdiff$Enzyme <- factor(SGdiff$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(SGdiff,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
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
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

####################################################################################################################
####################################################################################################################

#Prel enzymes: filter Prel enzymes from enz.genelist file and add the mean RPKM value at Prel
Enz_Prel <- enz.genelist[which(enz.genelist$Symbol %in% Prel_UP),]
Enz_Prel_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% Prel_UP,7:9])
Enz_Prel$RPKM <- Enz_Prel_RPKM[match(Enz_Prel$Symbol,names(Enz_Prel_RPKM))]

#construct table for dotplot for each mark
Enz_Prel_Hac <- filter(Enz_Prel,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_Prel$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_Prel_H3K4me <- filter(Enz_Prel,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_Prel$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_Prel_H3K36me <- filter(Enz_Prel,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_Prel$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_Prel_H3K9me <- filter(Enz_Prel,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_Prel$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_Prel_H3K27me <- filter(Enz_Prel,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_Prel$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_Prel_H4K20me <- filter(Enz_Prel,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_Prel$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_Prel_H3K79me <- filter(Enz_Prel,H3K79me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_Prel$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_Prel_H2AZK7me <- filter(Enz_Prel,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_Prel_H4K12me <- filter(Enz_Prel,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

Prel <- rbind(Enz_Prel_Hac,Enz_Prel_H3K4me,
                Enz_Prel_H3K36me, Enz_Prel_H3K9me,
                Enz_Prel_H3K27me, Enz_Prel_H4K20me,
                Enz_Prel_H3K79me, Enz_Prel_H2AZK7me,
                Enz_Prel_H4K12me) 

#create column with descending order to arrange graph
Prel$order <-  rep(nrow(Prel):1) 

#Factorize variables so they appear in the desired order
Prel$HK <-  factor(Prel$HK, 
                     levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

Prel$Enzyme <- factor(Prel$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(Prel,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(5,50,100,200,500),
                        limits = c(.05, 500),
                        labels = c("5","50","100","200","500"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

####################################################################################################################
####################################################################################################################

#LZ enzymes: filter LZ enzymes from enz.genelist file and add the mean RPKM value at LZ
Enz_LZ <- enz.genelist[which(enz.genelist$Symbol %in% LZ_UP),]
Enz_LZ_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% LZ_UP,10:12])
Enz_LZ$RPKM <- Enz_LZ_RPKM[match(Enz_LZ$Symbol,names(Enz_LZ_RPKM))]

#construct table for dotplot for each mark
Enz_LZ_Hac <- filter(Enz_LZ,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_LZ$Hac%in%ac_index))) %>%  rename("HPTM"= "Hac")

Enz_LZ_H3K4me <- filter(Enz_LZ,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_LZ$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_LZ_H3K36me <- filter(Enz_LZ,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_LZ$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_LZ_H3K9me <- filter(Enz_LZ,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_LZ$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_LZ_H3K27me <- filter(Enz_LZ,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_LZ$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_LZ_H4K20me <- filter(Enz_LZ,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_LZ$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_LZ_H3K79me <- filter(Enz_LZ,H3K79me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_LZ$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_LZ_H2AZK7me <- filter(Enz_LZ,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_LZ_H4K12me <- filter(Enz_LZ,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

LZ <- rbind(Enz_LZ_Hac,Enz_LZ_H3K4me,
              Enz_LZ_H3K36me, Enz_LZ_H3K9me,
              Enz_LZ_H3K27me, Enz_LZ_H4K20me,
              Enz_LZ_H3K79me, Enz_LZ_H2AZK7me,
              Enz_LZ_H4K12me) 

#create column with descending order to arrange graph
LZ$order <-  rep(nrow(LZ):1)

#Factorize variables so they appear in the desired order
LZ$HK <-  factor(LZ$HK, 
                   levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H4K12me","H3K79me","H3K9me","H3K27me","H4K20me"))

LZ$Enzyme <- factor(LZ$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(LZ,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(5,50,100,200,500),
                        limits = c(.05, 500),
                        labels = c("5","50","100","200","500"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

####################################################################################################################
####################################################################################################################

#PD enzymes: filter PD enzymes from enz.genelist file and add the mean RPKM value at PD
Enz_PD <- enz.genelist[which(enz.genelist$Symbol %in% PD_UP),]
Enz_PD_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% PD_UP,13:15])
Enz_PD$RPKM <- Enz_PD_RPKM[match(Enz_PD$Symbol,names(Enz_PD_RPKM))]

#construct table for dotplot for each mark
Enz_PD_Hac <- filter(Enz_PD,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_PD$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_PD_H3K4me <- filter(Enz_PD,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_PD$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_PD_H3K36me <- filter(Enz_PD,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_PD$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_PD_H3K9me <- filter(Enz_PD,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_PD$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_PD_H3K27me <- filter(Enz_PD,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_PD$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_PD_H4K20me <- filter(Enz_PD,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_PD$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_PD_H3K79me <- filter(Enz_PD,H3K79me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_PD$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_PD_H2AZK7me <- filter(Enz_PD,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_PD_H4K12me <- filter(Enz_PD,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

PD <- rbind(Enz_PD_Hac,Enz_PD_H3K4me,
            Enz_PD_H3K36me, Enz_PD_H3K9me,
            Enz_PD_H3K27me, Enz_PD_H4K20me,
            Enz_PD_H3K79me, Enz_PD_H2AZK7me,
            Enz_PD_H4K12me) 

#create column with descending order to arrange graph
PD$order <-  rep(nrow(PD):1) #create column with descending order to arrange graph

#Factorize variables so they appear in the desired order
PD$HK <-  factor(PD$HK, 
                 levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

PD$Enzyme <- factor(PD$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(PD,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(5,50,100,200,500),
                        limits = c(.05, 500),
                        labels = c("5","50","100","200","500"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

####################################################################################################################
####################################################################################################################

#RStid enzymes: filter RStid enzymes from enz.genelist file and add the mean RPKM value at RStid
Enz_RStid <- enz.genelist[which(enz.genelist$Symbol %in% RStid_UP),]
Enz_RStid_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% RStid_UP,16:18])
Enz_RStid$RPKM <- Enz_RStid_RPKM[match(Enz_RStid$Symbol,names(Enz_RStid_RPKM))]

#construct table for dotplot for each mark
Enz_RStid_Hac <- filter(Enz_RStid,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_RStid$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_RStid_H3K4me <- filter(Enz_RStid,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_RStid$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_RStid_H3K36me <- filter(Enz_RStid,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_RStid$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_RStid_H3K9me <- filter(Enz_RStid,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_RStid$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_RStid_H3K27me <- filter(Enz_RStid,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_RStid$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_RStid_H4K20me <- filter(Enz_RStid,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_RStid$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_RStid_H3K79me <- filter(Enz_RStid,H3K79me %in% me_index) %>%
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_RStid$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_RStid_H2AZK7me <- filter(Enz_RStid,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_RStid_H4K12me <- filter(Enz_RStid,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

RStid <- rbind(Enz_RStid_Hac,Enz_RStid_H3K4me,
            Enz_RStid_H3K36me, Enz_RStid_H3K9me,
            Enz_RStid_H3K27me, Enz_RStid_H4K20me,
            Enz_RStid_H3K79me,Enz_RStid_H2AZK7me,
            Enz_RStid_H4K12me) 

#create column with descending order to arrange graph
RStid$order <-  rep(nrow(RStid):1)

#Factorize variables so they appear in the desired order
RStid$HK <-  factor(RStid$HK, 
                 levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

RStid$Enzyme <- factor(RStid$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(RStid,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(5,50,100,200,400),
                        limits = c(.05, 400),
                        labels = c("5","50","100","200","400"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

####################################################################################################################
####################################################################################################################

#unnasigned enzymes: filter RStid enzymes from enz.genelist file and add the mean RPKM value at RStid
Enz_unnasigned <- enz.genelist[which(enz.genelist$Symbol %in% unassigned_enz),]
Enz_unnasigned_RPKM <- rowMeans(RPKMs[rownames(RPKMs) %in% unassigned_enz,1:15])
Enz_unnasigned$RPKM <- Enz_unnasigned_RPKM[match(Enz_unnasigned$Symbol,names(Enz_unnasigned_RPKM))]

#construct table for dotplot for each mark
Enz_unnasigned_Hac <- filter(Enz_unnasigned,Hac %in% ac_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,Hac)) %>% mutate(HK=rep("Hac",sum(Enz_unnasigned$Hac%in%ac_index))) %>% rename("HPTM"= "Hac")

Enz_unnasigned_H3K4me <- filter(Enz_unnasigned,H3K4me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K4me)) %>% mutate(HK=rep("H3K4me",sum(Enz_unnasigned$H3K4me%in%me_index))) %>% rename("HPTM"= "H3K4me")

Enz_unnasigned_H3K36me <- filter(Enz_unnasigned,H3K36me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K36me)) %>% mutate(HK=rep("H3K36me",sum(Enz_unnasigned$H3K36me%in%me_index))) %>% rename("HPTM"= "H3K36me")

Enz_unnasigned_H3K9me <- filter(Enz_unnasigned,H3K9me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K9me)) %>% mutate(HK=rep("H3K9me",sum(Enz_unnasigned$H3K9me%in%me_index))) %>% rename("HPTM"= "H3K9me")

Enz_unnasigned_H3K27me <- filter(Enz_unnasigned,H3K27me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K27me)) %>% mutate(HK=rep("H3K27me",sum(Enz_unnasigned$H3K27me%in%me_index))) %>% rename("HPTM"= "H3K27me")

Enz_unnasigned_H4K20me <- filter(Enz_unnasigned,H4K20me %in% me_index) %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K20me)) %>% mutate(HK=rep("H4K20me",sum(Enz_unnasigned$H4K20me%in%me_index))) %>% rename("HPTM"= "H4K20me")

Enz_unnasigned_H3K79me <- filter(Enz_unnasigned,H3K79me %in% me_index) %>%
  dplyr::select(c(RPKM,Symbol,Enzyme,H3K79me)) %>% mutate(HK=rep("H3K79me",sum(Enz_unnasigned$H3K79me%in%me_index))) %>% rename("HPTM"= "H3K79me")

Enz_unnasigned_H2AZK7me <- filter(Enz_unnasigned,H2AZK7me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H2AZK7me)) %>% mutate(HK=rep("H2AZK7me",1)) %>% rename("HPTM"= "H2AZK7me")

Enz_unnasigned_H4K12me <- filter(Enz_unnasigned,H4K12me=="me1") %>% 
  dplyr::select(c(RPKM,Symbol,Enzyme,H4K12me)) %>% mutate(HK=rep("H4K12me",1)) %>% rename("HPTM"= "H4K12me")

unnasigned <- rbind(Enz_unnasigned_Hac,Enz_unnasigned_H3K4me,
               Enz_unnasigned_H3K36me, Enz_unnasigned_H3K9me,
               Enz_unnasigned_H3K27me, Enz_unnasigned_H4K20me,
               Enz_unnasigned_H3K79me,Enz_unnasigned_H2AZK7me,
               Enz_unnasigned_H4K12me) 

#create column with descending order to arrange graph
unnasigned$order <-  rep(nrow(unnasigned):1)

#Factorize variables so they appear in the desired order
unnasigned$HK <-  factor(unnasigned$HK, 
                    levels = c("Hac","H3K4me","H3K36me","H2AZK7me","H3K79me","H3K9me","H3K27me","H4K20me"))

unnasigned$Enzyme <- factor(unnasigned$Enzyme, levels = c("Writer","Eraser"))

#create plot     
ggplot(unnasigned,
       aes(x=HK,y=reorder(Symbol,order),
           color=Enzyme))+
  geom_point(aes(size=RPKM))+
  scale_color_manual(name="Enzyme",
                     labels=c("Writer","Eraser"),
                     values = c("blue","red"))+
  scale_size_continuous(name = "RPKM",
                        breaks = c(1,10,50,100,200),
                        limits = c(.05, 400),
                        labels = c("1","10","50","100","200"),
                        range=c(1,10))+
  theme_gray()+
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  xlab("")+ylab("")+
  geom_text_repel(aes(label=HPTM), size=2.5, nudge_x=0.5)

