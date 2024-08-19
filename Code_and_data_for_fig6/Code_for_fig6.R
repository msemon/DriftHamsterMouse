library(ggplot2)
library(cowplot)

id="Shh"

# For panel A
sp_names <- c(`bat` = "bat",
                `mus` = "mouse")
to_names <- c(`fore.limb` = "fore.limb",
                `hind.limb` = "hind.limb")

models <-read.table(file=paste0(id,"_2_species_data_for_Fig6A.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig6A.txt"),header=T,sep="\t")
 
p=ggplot(dats, aes(x=times, y=value, group=species)) + 
facet_grid(species ~ tissue, labeller = labeller(species=sp_names,tissue=to_names))+theme_bw() +
geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
scale_colour_manual(name="Model",
                        values=c("Simple"="grey",
                                 "bat-mouse"="blue",
                                 "hind-fore"="red",
                                 "4 curves"="black"))+
ggtitle(paste0(id))+
    theme(plot.title = element_text(size = 12, face = "bold"))

save_plot("fig6_panelA.pdf", ncol=2,nrow=3,
            p ,limitsize = FALSE )     
  

# For panel B
id="Shh"


models <-read.table(file=paste0(id,"_2_species_data_for_Fig6B.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig6B.txt"),header=T,sep="\t")

sp_names <- c(`bat` = "bat",
              `mus` = "mouse")
to_names <- c(`fore-limb` = "fore.limb",
              `hind-limb` = "hind.limb")

p=ggplot(dats, aes(x=times, y=value, group=species,col=tissue)) + 
  facet_grid(species ~ ., labeller = labeller(species=sp_names,tissue=to_names))+theme_bw() +
  geom_point(size=4,alpha = 0.5) + theme_bw(base_size = 22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
  geom_line(data=models,aes(x=times, y=value,col=tissue,group = tissue))+
  scale_colour_manual(name="Limb",
                      values=c("fore-limb"="black",
                               "hind-limb"="grey"))+
  ggtitle(id)+
  theme(plot.title = element_text(size = 12, face = "bold"))
save_plot("fig6_panelB.pdf", ncol=1.5,nrow=2,
          p ,limitsize = FALSE )    




# For panel C
id="Fgf8"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig6C.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig6C.txt"),header=T,sep="\t")

sp_names <- c(`bat` = "bat",
              `mus` = "mouse")
to_names <- c(`fore-limb` = "fore.limb",
              `hind-limb` = "hind.limb")

p=ggplot(dats, aes(x=times, y=value, group=species,col=tissue)) + 
  facet_grid(species ~ ., labeller = labeller(species=sp_names,tissue=to_names))+theme_bw() +
  geom_point(size=4,alpha = 0.5) + theme_bw(base_size = 22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
  geom_line(data=models,aes(x=times, y=value,col=tissue,group = tissue))+
  scale_colour_manual(name="Limb",
                      values=c("fore-limb"="black",
                               "hind-limb"="grey"))+
  ggtitle(id)+
  theme(plot.title = element_text(size = 12, face = "bold"))
save_plot("fig6_panelC.pdf", ncol=1.5,nrow=2,
          p ,limitsize = FALSE )    



# For panel D
id="Grem1"

models <-read.table(file=paste0(id,"_2_species_data_for_Fig6D.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_data_for_Fig6D.txt"),header=T,sep="\t")

sp_names <- c(`bat` = "bat",
              `mus` = "mouse")
to_names <- c(`fore-limb` = "fore.limb",
              `hind-limb` = "hind.limb")

p=ggplot(dats, aes(x=times, y=value, group=species,col=tissue)) + 
  facet_grid(species ~ ., labeller = labeller(species=sp_names,tissue=to_names))+theme_bw() +
  geom_point(size=4,alpha = 0.5) + theme_bw(base_size = 22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
  geom_line(data=models,aes(x=times, y=value,col=tissue,group = tissue))+
  scale_colour_manual(name="Limb",
                      values=c("fore-limb"="black",
                               "hind-limb"="grey"))+
  ggtitle(id)+
  theme(plot.title = element_text(size = 12, face = "bold"))
save_plot("fig6_panelD.pdf", ncol=1.5,nrow=2,
          p ,limitsize = FALSE )     
