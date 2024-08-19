
library(ggplot2)
library(viridis)
library(cowplot)




colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")


sp_names <- c(`ham` = "hamster",
              `mus` = "mouse")
to_names <- c(`md` = "lower",
              `mx` = "upper")
id="Osr2"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig.txt"),header=T,sep="\t")

p<-ggplot(dats, aes(x=times, y=value, group= paste(species,jaw), col = paste(species,jaw))) + 
  facet_grid(species ~., labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
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

ggsave("figS3_panelA.pdf",p,width = 12, height = 10, units = "cm") 
  

id="Sfrp2"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig.txt"),header=T,sep="\t")

p<-ggplot(dats, aes(x=times, y=value, group= paste(species,jaw), col = paste(species,jaw))) + 
  facet_grid(species ~., labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
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

ggsave("figS3_panelB.pdf",p,width = 12, height = 10, units = "cm") 


## Bottom
dfgene<-read.table("data_BL_gene_Sfrp2.txt",h=T)
g=ggplot(dfgene, aes(y=tooth, fill=gene, x= tis)) + geom_tile(color="white", size=0.1)+
    scale_y_discrete(name ="Tooth",labels=c("md"="lower","mx"="upper"))+
    scale_x_discrete(name ="Tooth Part",labels=c("B"="buccal","L"="lingual"))+
    scale_fill_viridis(name=id,direction=-1)+
    coord_equal()+
    theme(axis.ticks=element_blank())+
    theme(legend.text=element_text(size=6))+ ggpubr::rotate_x_text()
  ggsave("FigS3_panel_bottom_Sfrp2.pdf",g)

  
  dfgene<-read.table("data_BL_gene_Osr2.txt",h=T)
  g=ggplot(dfgene, aes(y=tooth, fill=gene, x= tis)) + geom_tile(color="white", size=0.1)+
    scale_y_discrete(name ="Tooth",labels=c("md"="lower","mx"="upper"))+
    scale_x_discrete(name ="Tooth Part",labels=c("B"="buccal","L"="lingual"))+
    scale_fill_viridis(name=id,direction=-1)+
    coord_equal()+
    theme(axis.ticks=element_blank())+
    theme(legend.text=element_text(size=6))+ ggpubr::rotate_x_text()
  ggsave("FigS3_panel_bottom_Osr2.pdf",g)
  
  