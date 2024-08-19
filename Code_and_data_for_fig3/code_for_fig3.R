library(ggplot2)
library(ggpubr)
library(viridis)
library(cowplot)

#### For Panel A

# temporary: for converting lightseet dev time to relative time
#MusZero=14.65549
#HamZero=12.30508
#MusTen=17.98604
#HamTen=14.57598
#10*(15.1938006194863-MusZero)/(MusTen-MusZero)
#10*(12.6375183963396-HamZero)/(HamTen-HamZero)

  
tabmes<-read.table(file="tab_mes_for_fig3_panelA.txt",sep="\t",h=T)
tabmes$type=ifelse(grepl("lightsheet",tabmes$sample),"lightsheet","deconvolution")

p=ggplot(tabmes, aes(x=timerelGAM, y=median*100, color=tooth)) + 
  facet_grid(species~type)+theme_bw()+
  geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
  ylim(0,100)+
  geom_segment(data=tabmes,aes(x=timerelGAM, xend=timerelGAM,y=mini*100,yend=maxi*100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time",y="% mesenchyme")+
  scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename="fig3_panelA.pdf",ncol=2,nrow=2,p)


#### For Panel C

tabRoma<-read.table(file="tab_roma_for_fig3_panelC.txt",sep="\t",h=T)


g<-ggplot(tabRoma, aes(y=tooth, fill=mean_activ, x= type)) + geom_tile(color="white", size=0.1)+
  scale_y_discrete(name ="Tooth",labels=c("md"="lower","mx"="upper"))+
  scale_x_discrete(name ="Tooth Part",labels=c("B"="buccal","L"="lingual"))+
  scale_fill_viridis(name="Activation",direction=-1)+
  coord_equal()+
  facet_grid( ~ pathComp) +
  theme(axis.ticks=element_blank())+
  theme(legend.text=element_text(size=6))+ ggpubr::rotate_x_text()

ggsave("fig3_panelC.pdf",g,width=20,height=10,unit="cm")


#### For Panel D

tablingual<-read.table(file="tab_lingual_for_fig3_panelD.txt",sep="\t",h=T)

p=ggplot(tablingual, aes(x=timerelGAM, y=Lmedian*100, color=tooth)) + 
  facet_grid(. ~ species, scales = "free_x")+theme_bw()+
  geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
  ylim(0,100)+
  geom_segment(data=tablingual,aes(x=timerelGAM, xend=timerelGAM,y=Lmini*100,yend=Lmaxi*100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time",y="% lingual")+
  scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename="fig3_panelD.pdf",ncol=2,nrow=1,p)