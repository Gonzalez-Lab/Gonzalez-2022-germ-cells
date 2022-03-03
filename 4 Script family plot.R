library(tidyverse)

load("Enzymes MKs RPKM matrix.Rdata")

# define histones PTMs indexes
ac_index <- c("Hac","H3ac","H4ac","H2/H4ac","H3/H4ac") # 1=H2/3/4, 3=H3, 4=H4, 24= H2/4, 34 = H3/4
me_index <- c("me1","me2","me3") #me1 me2 me3, indicates highest position described for enzyme


## Tables and plots for histone acetylation Hac

Hac_families <- Enz_MK_RPKM %>% 
  filter(Hac %in% ac_index) %>% 
  group_by(Family) %>% summarize(sumSGThy1=sum(SGThy1),
                                 sumSGKit=sum(SGKit),
                                 sumSCyte=sum(SCyte),
                                 sumSTid=sum(STid),
                                 sumSperm=sum(Sperm))

Hac_families<- as.data.frame(Hac_families)
Hac_families$Enzyme <- c(rep("Eraser",4),rep("Writer",6))

#Hac eraser N
n=6
Hac_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                          rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                   Family=rep(Hac_families[5:10,1],5),
                                   RPKM= c(Hac_families[5:10,2],
                                           Hac_families[5:10,3],
                                           Hac_families[5:10,4],
                                           Hac_families[5:10,5],
                                           Hac_families[5:10,6])))



Hac_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                             y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
 theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "Hac writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,800))+
  xlab("")+ylab("")
                   
#Hac eraser N
n=4

Hac_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                          rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                   Family=rep(Hac_families[1:4,1],5),
                                   RPKM= c(Hac_families[1:4,2],
                                           Hac_families[1:4,3],
                                           Hac_families[1:4,4],
                                           Hac_families[1:4,5],
                                           Hac_families[1:4,6])))



Hac_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                             y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "Hac erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,800))+
  xlab("")+ylab("")

#########################################################################################
## Tables and plots for H3K4me1/2/3 methylation

H3K4me_families <- Enz_MK_RPKM %>% 
  filter(H3K4me %in% me_index[-1]) %>% mutate(Family_Mark=paste(Family,H3K4me)) %>% 
  group_by(Family_Mark) %>% summarize(sumSGThy1=sum(SGThy1),
                                 sumSGKit=sum(SGKit),
                                 sumSCyte=sum(SCyte),
                                 sumSTid=sum(STid),
                                 sumSperm=sum(Sperm))

H3K4me_families<- as.data.frame(H3K4me_families)
H3K4me_families$Enzyme <- c(rep("Eraser",3),rep("Writer",5))

#H3K4me writers N
n=5
H3K4me_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                          rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                   Family=rep(H3K4me_families[4:8,1],5),
                                   RPKM= c(H3K4me_families[4:8,2],
                                           H3K4me_families[4:8,3],
                                           H3K4me_families[4:8,4],
                                           H3K4me_families[4:8,5],
                                           H3K4me_families[4:8,6])))



H3K4me_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                             y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K4me2/3 writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,300))+
  xlab("")+ylab("")

#H3K4me eraser N
n=3

H3K4me_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                          rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                   Family=rep(H3K4me_families[1:3,1],5),
                                   RPKM= c(H3K4me_families[1:3,2],
                                           H3K4me_families[1:3,3],
                                           H3K4me_families[1:3,4],
                                           H3K4me_families[1:3,5],
                                           H3K4me_families[1:3,6])))



H3K4me_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                             y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K4me2/3 erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,300))+
  xlab("")+ylab("")

#########################################################################################
## Tables and plots for H3K36me2/3 methylation

H3K36me_families <- Enz_MK_RPKM %>% 
  filter(H3K36me %in% me_index) %>% mutate(Family_Mark=paste(Family,H3K36me)) %>% 
  group_by(Family_Mark) %>% summarize(sumSGThy1=sum(SGThy1),
                                 sumSGKit=sum(SGKit),
                                 sumSCyte=sum(SCyte),
                                 sumSTid=sum(STid),
                                 sumSperm=sum(Sperm))

H3K36me_families<- as.data.frame(H3K36me_families)
H3K36me_families$Enzyme <- c(rep("Eraser",4),rep("Writer",6))

#H3K36me writers N
n=6
H3K36me_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                             rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                      Family=rep(H3K36me_families[5:10,1],5),
                                      RPKM= c(H3K36me_families[5:10,2],
                                              H3K36me_families[5:10,3],
                                              H3K36me_families[5:10,4],
                                              H3K36me_families[5:10,5],
                                              H3K36me_families[5:10,6])))

H3K36me_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K36me2/3 writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,320))+
  xlab("")+ylab("")

#H3K36me eraser N
n=4

H3K36me_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                             rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                      Family=rep(H3K36me_families[1:4,1],5),
                                      RPKM= c(H3K36me_families[1:4,2],
                                              H3K36me_families[1:4,3],
                                              H3K36me_families[1:4,4],
                                              H3K36me_families[1:4,5],
                                              H3K36me_families[1:4,6])))



H3K36me_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K36me2/3 erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,320))+
  xlab("")+ylab("")

#########################################################################################
## Tables and plots for H3K9me2/3 methylation

H3K9me_families <- Enz_MK_RPKM %>% 
  filter(H3K9me %in% me_index[-1]) %>% mutate(Family_Mark=paste(Family,H3K9me)) %>% 
  group_by(Family_Mark) %>% summarize(sumSGThy1=sum(SGThy1),
                                      sumSGKit=sum(SGKit),
                                      sumSCyte=sum(SCyte),
                                      sumSTid=sum(STid),
                                      sumSperm=sum(Sperm))

H3K9me_families<- as.data.frame(H3K9me_families)
H3K9me_families$Enzyme <- c(rep("Eraser",4),rep("Writer",3))

#H3K9me writers N
n=3
H3K9me_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                              rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                       Family=rep(H3K9me_families[5:7,1],5),
                                       RPKM= c(H3K9me_families[5:7,2],
                                               H3K9me_families[5:7,3],
                                               H3K9me_families[5:7,4],
                                               H3K9me_families[5:7,5],
                                               H3K9me_families[5:7,6])))



H3K9me_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                 y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K9me2/3 writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,460))+
  xlab("")+ylab("")

#H3K9me eraser N
n=4
H3K9me_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                              rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                       Family=rep(H3K9me_families[1:4,1],5),
                                       RPKM= c(H3K9me_families[1:4,2],
                                               H3K9me_families[1:4,3],
                                               H3K9me_families[1:4,4],
                                               H3K9me_families[1:4,5],
                                               H3K9me_families[1:4,6])))



H3K9me_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                 y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K9me2/3 erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,460))+
  xlab("")+ylab("")

###########################################################################################################
## Tables and plots for H3K27me2/3 methylation

H3K27me_families <- Enz_MK_RPKM %>% 
  filter(H3K27me %in% me_index) %>% mutate(Family_Mark=paste(Family,H3K27me)) %>% 
  group_by(Family_Mark) %>% summarize(sumSGThy1=sum(SGThy1),
                                      sumSGKit=sum(SGKit),
                                      sumSCyte=sum(SCyte),
                                      sumSTid=sum(STid),
                                      sumSperm=sum(Sperm))

H3K27me_families<- as.data.frame(H3K27me_families)
H3K27me_families$Enzyme <- c(rep("Eraser",2),rep("Writer",1))

#H3K27me writers N
n=1
H3K27me_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                             rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                      Family=rep(H3K27me_families[3,1],5),
                                      RPKM= c(H3K27me_families[3,2],
                                              H3K27me_families[3,3],
                                              H3K27me_families[3,4],
                                              H3K27me_families[3,5],
                                              H3K27me_families[3,6])))



H3K27me_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K27me2/3 writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,120))+
  xlab("")+ylab("")

#H3K27me eraser N
n=2

H3K27me_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                             rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                      Family=rep(H3K27me_families[1:2,1],5),
                                      RPKM= c(H3K27me_families[1:2,2],
                                              H3K27me_families[1:2,3],
                                              H3K27me_families[1:2,4],
                                              H3K27me_families[1:2,5],
                                              H3K27me_families[1:2,6])))



H3K27me_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H3K27me2/3 erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,125))+
  xlab("")+ylab("")

###########################################################################################################
## Tables and plots for H4K20me2/3 methylation

H4K20me_families <- Enz_MK_RPKM %>% 
  filter(H4K20me %in% me_index[-1]) %>% mutate(Family_Mark=paste(Family,H4K20me)) %>% 
  group_by(Family_Mark) %>% summarize(sumSGThy1=sum(SGThy1),
                                      sumSGKit=sum(SGKit),
                                      sumSCyte=sum(SCyte),
                                      sumSTid=sum(STid),
                                      sumSperm=sum(Sperm))

H4K20me_families<- as.data.frame(H4K20me_families)
H4K20me_families$Enzyme <- c(rep("Eraser",2),rep("Writer",4))

#H4K20me writers N
n=4
H4K20me_writers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                              rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                       Family=rep(H4K20me_families[3:6,1],5),
                                       RPKM= c(H4K20me_families[3:6,2],
                                               H4K20me_families[3:6,3],
                                               H4K20me_families[3:6,4],
                                               H4K20me_families[3:6,5],
                                               H4K20me_families[3:6,6])))



H4K20me_writers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                 y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H4K20me writers")+
  scale_fill_brewer(palette="Blues", direction = -1)+
  ylim(c(0,120))+
  xlab("")+ylab("")

#H4K20me eraser N
n=2

H4K20me_erasers <- as.data.frame(cbind(Cell=c(rep("SGThy1",n),rep("SGKit",n),
                                              rep("SCyte",n),rep("STid",n),rep("Sperm",n)),
                                       Family=rep(H4K20me_families[1:2,1],5),
                                       RPKM= c(H4K20me_families[1:2,2],
                                               H4K20me_families[1:2,3],
                                               H4K20me_families[1:2,4],
                                               H4K20me_families[1:2,5],
                                               H4K20me_families[1:2,6])))



H4K20me_erasers  %>%  ggplot(aes(x=factor(Cell, levels = c("SGThy1","SGKit","SCyte","STid","Sperm")), 
                                 y=as.numeric(RPKM), fill=Family))+
  geom_bar(stat="identity", position = "stack")+
  theme_classic()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  labs(title = "H4K20me erasers")+
  scale_fill_brewer(palette="Reds", direction = -1)+
  ylim(c(0,125))+
  xlab("")+ylab("")
