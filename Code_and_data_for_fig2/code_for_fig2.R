library(ggplot2)
library(ggrepel)

# Panel A


theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


datacusps <- read.table("datasamples_for_fig2A.txt",h=T)
g=ggplot(datacusps, aes(y= factor(species), x=time, 
                        label=as.character(morphostage), 
                        color=type))+
  xlab("time")+
  ylab("species")+
  geom_point(position=position_dodge(.5), size=4, shape=20) +
  geom_errorbar(aes(xmin=timeL, xmax=timeU),width=.01,
                position=position_dodge(.5)) + 
  geom_label_repel(force=0.1,force_pull=1,direction="y")+
  #coord_flip()+
  theme_dotplot+
  scale_color_manual(values=c("#E69F00", "#56B4E9","#999999"))  
ggsave(g,file="fig2_panelA.pdf",width = 15, height = 15, units = "cm")


# Panel B


colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")

c=read.table("data_1_for_pca_fig2B.txt",h=T)
eig=read.table("data_2_for_pca_fig2B.txt",h=T)$x
eig=round(eig,1)
# here it is possible to change i and j to see other PCA axes
i=2
j=1
n1=paste("Axis",i,sep="")
n2=paste("Axis",j,sep="")

g=ggplot(c, aes(x=c[,n1], y=c[,n2], color=as.numeric(stagenorm), shape = paste(tooth,species,sep="."))) +geom_point(size=5)+
  #ggplot(c, aes(x=c[,n1], y=c[,n2], color=paste(espece,stade), shape = machoire)) +geom_point(size=1) + 
  theme_bw(base_size = 12)+
  ylim(-100,100)+
  xlim(-100,100)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  labs(x=paste(n1,eig[i],"%"),y=paste(n2,eig[j],"%"),shape="tooth.sp", colour="Dev. time")+
  scale_shape_manual(name="Tooth",
                     values=c("md.ham"=16,
                              "mx.ham"=15,
                              "md.mus"=17,
                              "mx.mus"=18),
                     labels=c("md.ham"= "lower hamster",
                              "mx.ham"="upper hamster",
                              "md.mus"="lower mus",
                              "mx.mus"="upper mus"))
ggsave(g,file="fig2_panelB.pdf",width = 15, height = 15, units = "cm")


# Panel C

id="Fgf10"

models <-read.table(file=paste0(id,"_2_species_for_Fig2C.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_for_Fig2C.txt"),header=T,sep="\t")

sp_names <- c(`ham` = "hamster",
              `mus` = "mouse")
to_names <- c(`md` = "lower",
              `mx` = "upper")


colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")

p<-ggplot(dats, aes(x=times, y=value, group= paste(species,jaw), col = paste(species,jaw))) + 
  facet_grid(species ~jaw, labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
  geom_point(size=2,alpha = 0.4) + theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
  geom_line(data=models,aes(x=times, y=value,group= paste(species,jaw),col = paste(species,jaw)))+
  scale_colour_manual(name="tooth",
                      values=colo,
                      labels=c("ham md"="ham lower",
                               "ham mx"="ham upper",
                               "mus md"="mus lower",
                               "mus mx"="mus upper"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), strip.background = element_blank(), strip.text.x = element_blank())#+

ggsave("fig2_panelC.pdf",p,width = 15, height = 10, units = "cm") 


# Panel D

metaData_dist<-read.table(file="data_distances_for_fig2D.txt",h=T,sep="\t")


theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())

w10<-ggplot(metaData_dist[metaData_dist$typedist=="dist",],aes(x=window,y=mD,shape=comparison))+
  geom_point()+theme_dotplot+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  scale_shape_manual(values=c(1, 10, 4,3))

ggsave(w10,file="fig2_panelD.pdf",w = 14, h = 15, units="cm")
