library(ggplot2)
library(dplyr)
library(enrichplot)

## Panel A

gorich_BP_data=readRDS("coevo_BP_table.rds")
edo <- pairwise_termsim(gorich_BP_data)
p<-emapplot(edo, layout.params = list(layout="kk"), cluster.params =list(cluster = T, n=10, legend = T), cex_label_group = 1.5, node_label = "group")
ggsave("FigSCoevo_panelA_v1.pdf",p,width=15,height=15)



gopatp = gorich_BP_data  %>% 
  mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio)) %>% 
  filter(p.adjust < 0.1)  %>%
  ggplot(showCategory = 100,
         aes(FoldEnrichment, forcats::fct_reorder(Description, FoldEnrichment ))) + 
  geom_segment(aes(xend=1, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("FoldEnrichment") +
  ylab(NULL) + 
  ggtitle("Biological Process enrichment")
ggsave(paste0("FigSCoevo_panelA_v2.pdf"),gopatp,width=15,height=20)



## Panel B

gopat=read.table("coevo_pathwayenrichment_table_FE.txt",sep="\t")


gopatp<-gopat%>%
  filter(p.adjust < 0.1)  %>%
  ggplot(showCategory = 100,
         aes(FoldEnrichment, forcats::fct_reorder(Description, FoldEnrichment ))) + 
  geom_segment(aes(xend=1, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("FoldEnrichment") +
  ylab(NULL) + 
  ggtitle("pathway enrichment")
ggsave(paste0("FigSCoevo_panelB.pdf"),gopatp,width=15,height=20)


## Panel C


colo=c("ham md"="dark grey",
       "ham mx"="black",
       "mus md"="dark grey",
       "mus mx"="black")


sp_names <- c(`ham` = "hamster",
              `mus` = "mouse")
to_names <- c(`md` = "lower",
              `mx` = "upper")



id="Igf1"

models <-read.table(file=paste0(id,"_2_species_for_Fig.txt"),header=T,sep="\t")
dats <-read.table(file=paste0(id,"_1_species_for_Fig.txt"),header=T,sep="\t")

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

ggsave("FigSCoevo_panelC.pdf",p,width = 12, height = 10, units = "cm") 

