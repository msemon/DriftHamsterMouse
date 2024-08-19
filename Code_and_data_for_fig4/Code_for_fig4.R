library(ggplot2)
library(cowplot)
library(viridis)

### Panel A

# exemple (left)

id="Bmp4"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig4A.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig4A.txt"),header=T,sep="\t")

sp_names <- c(`ham` = "hamster",
              `mus` = "mouse")
to_names <- c(`md` = "lower",
              `mx` = "upper")
p=ggplot(dats, aes(x=times, y=value, group=species)) + 
  facet_grid(species ~ jaw, labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
  geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
  geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
  scale_colour_manual(name="Model",
                      values=c("Simple"="grey",
                               "Hamster-mouse"="blue",
                               "Upper-lower"="red",
                               "4 teeth"="black"))+
  ggtitle(id)+
  theme(plot.title = element_text(size = 12, face = "bold"))


save_plot("fig4_panelA_left.pdf", ncol=2,nrow=3,
          p ,limitsize = FALSE )   


# pathway coevolution (right)

tabcoevo <- read.table(file="tab_coevo_pathway_for_Fig4A.txt",sep="\t",h=T)


theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())

propclass3=ggplot(tabcoevo, aes(x = subset2, y=percent)) +
  geom_bar(stat="identity", width=0.5)+
  scale_x_discrete(labels=tabcoevo$subset3)+
  theme_dotplot +
  xlab("")+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle = 45))+
  coord_flip()

propclass3 = propclass3 + theme(legend.position = "none")+
  annotate("rect", ymin = 0, ymax = max(tabcoevo$percent)+0.1,  xmax = 13.2, xmin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(tabcoevo$percent)+0.1,  xmax = 11.2, xmin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(tabcoevo$percent)+0.1,  xmax = 9.2, xmin = 4.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(tabcoevo$percent)+0.1,  xmax = 9.2, xmin = 4.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(tabcoevo$percent)+0.1,  xmax = 4.2, xmin = 0.8, fill = "lightpink", alpha = 0.2)

ggsave("fig4_panelA_right.pdf",propclass3,width = 12, height = 10, units = "cm") 



# Panel B

id="Bmp4"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig4B.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig4B.txt"),header=T,sep="\t")


colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")

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

ggsave("fig4_panelB.pdf",p,width = 12, height = 10, units = "cm") 
  



# Panel C

id="Wif1"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig4C.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig4C.txt"),header=T,sep="\t")

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

ggsave("fig4_panelC.pdf",p,width = 12, height = 10, units = "cm") 




# Panel D


id="Dlx1"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig4D.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig4D.txt"),header=T,sep="\t")

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

ggsave("fig4_panelD.pdf",p,width = 12, height = 10, units = "cm") 



# Panel E

## left


theme_barplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())



dats <-read.table(file="datevo_for_fig4E1.txt",header=T,sep="\t")
dats$subset3=factor(dats$subset3,levels=c("eda (8)","notch (9)","tgfb (11)",
                                          "activin (16)","HH (29)","FGF (63)","Bmp (65)","Wnt (75)",
                                          "Pathways (247)","Keystone (79)","Dispensable (88)",
                                          "bite-it (213)","Total (13278)"))

propclass=ggplot(dats, aes(x = freq, y = subset3,fill = type)) + 
  theme_barplot +
  geom_bar(width = 0.5,position=position_dodge(),stat="identity")+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_manual(name="Model", values = c("Lower" = "purple", "Upper" = "darkgreen"),labels = c("Lower 2curves","Upper 2curves"))
#scale_fill_viridis(name="Model",discrete = TRUE, option = "D",labels = c("Lower 2curves","Upper 2curves"))

ggsave("fig4_panelE_left.pdf",propclass,width = 20, height = 10, units = "cm") 



## right

dats <-read.table(file="datevo_for_fig4E2.txt",header=T,sep="\t")

dats$subset3=factor(dats$subset3,levels=c("eda (8)","notch (9)","tgfb (11)",
                                          "activin (16)","HH (29)","FGF (63)","Bmp (65)","Wnt (75)",
                                          "Pathways (247)","Keystone (79)","Dispensable (88)",
                                          "bite-it (213)","Total (13278)"))


propclass=ggplot(dats, aes(x = freq, y = subset3,fill = variable)) +
  geom_bar(stat="identity")+ theme_barplot +
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "D")

ggsave("fig4_panelE_right.pdf",propclass,width = 20, height = 10, units = "cm") 





