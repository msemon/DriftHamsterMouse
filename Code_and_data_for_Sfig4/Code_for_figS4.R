
library(ggplot2)
library(cowplot)




colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")


sp_names <- c(`ham` = "hamster",
              `mus` = "mouse")
to_names <- c(`md` = "lower",
              `mx` = "upper")


# Panel A


id="Bmp4"

models <-read.table(file=paste0(id,"_2_species_for_FigS4.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_for_FigS4.txt"),header=T,sep="\t")

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

ggsave("figS4_panelA.pdf",p,width = 12, height = 10, units = "cm") 




# Panel B

id="Dkk1"

models <-read.table(file=paste0(id,"_2_species_for_FigS4.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_for_FigS4.txt"),header=T,sep="\t")

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

ggsave("figS4_panelB.pdf",p,width = 12, height = 10, units = "cm") 




# Panel C

id="Wif1"

models <-read.table(file=paste0(id,"_2_species_for_FigS4.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_for_FigS4.txt"),header=T,sep="\t")

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

ggsave("figS4_panelC.pdf",p,width = 12, height = 10, units = "cm") 

