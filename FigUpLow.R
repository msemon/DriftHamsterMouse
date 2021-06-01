##### FIGURE UPPER/LOWER


setwd("~/Documents/Projet_Drift/MappingGenome/Counts/Spline_Final/TimeGAM")
options(connectionObserver = NULL)
library(org.Mm.eg.db)
library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(ade4)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(forcats)
library(gprofiler2)
library(DOSE)
library(clusterProfiler)
library(GOSemSim)
library(viridis)
library(readxl)
library(enrichplot)


resmod3=data.frame(mus.padj=res_mus_vs_tooth[[2]]$padj,
                  # mus.diste=res_mus_vs_tooth[[2]]$distearly,
                  # mus.distm=res_mus_vs_tooth[[2]]$distmid,
                  # mus.distl=res_mus_vs_tooth[[2]]$distlate,
                   ham.padj=res_ham_vs_tooth[[2]]$padj)
                 #  ham.diste=res_ham_vs_tooth[[2]]$distearly,
                #   ham.distm=res_ham_vs_tooth[[2]]$distmid,
                #   ham.distl=res_ham_vs_tooth[[2]]$distlate)
row.names(resmod3)=row.names(res_mus_vs_tooth[[2]])


resmod3$mus="1 curve"
resmod3$mus[resmod3$mus.padj<0.05]="Tooth"
resmod3$ham="1 curve"
resmod3$ham[resmod3$ham.padj<0.05]="Tooth"

table(mus=resmod3$mus,ham=resmod3$ham)

### 

pathway = read_excel("../pathways-Margaux-corMS.xlsx",sheet=1)
pathway=pathway[,1:5]

pathwaylist=unique(pathway$Symbol)
pathwaylistspe=lapply(split(pathway,pathway$Pathway),function(x){unique(x$Symbol)})
pathwaylist[pathwaylist%in%row.names(resnested)]


resmod3$model="1 curve"
resmod3$model[resmod3$mus=="Tooth"]="mus 2 teeth"
resmod3$model[resmod3$ham=="Tooth"]="ham 2 teeth"
resmod3$model[(resmod3$ham=="Tooth")&(resmod3$mus=="Tooth")]="2 teeth"

 
 getclassesmod3=function(list=biteit$V1,suf="bite-it",threshold=0.05)
 {
   subset=resmod3[row.names(resmod3)%in%list,]
   out=c(table(subset$model),nrow(subset),suf)
   names(out)=c(names(table(subset$model)),"N","subset")
   out=out[c("subset","N","1 curve","mus 2 teeth","ham 2 teeth","2 teeth")]
   return(out)
 }
 

n1=getclassesmod3(list=biteit$V1,suf="bite-it")
n2=getclassesmod3(list=dispensable$Gene..name,suf="Dispensable")
n3=getclassesmod3(list=keystone$Gene.name,suf="Keystone")
n6=getclassesmod3(list=row.names(resmod3),suf="Total")
n7=getclassesmod3(list=pathwaylist,suf="Pathways")


datvar=data.frame(rbind(n1,n2,n3,n6,n7))
#datvar=data.frame(rbind(n1,n5,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
datvarmelt$N=as.numeric(as.character(datvarmelt$N))

datvarmelt$subset=factor(datvarmelt$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))

datvarmelt2 <- datvarmelt %>%
  group_by(subset) %>%
  mutate(percent = value / sum(value)) %>%
  ungroup() %>%
  mutate(subset2 = paste0(subset," (",N,")"))%>%
  mutate(subset2 = fct_reorder(subset2, as.numeric(as.character(N))))


lev= levels(datvarmelt2$subset2)
lev1=lev
lev1[1]=lev[4]
lev1[4]=lev[3]
lev1[2:3]=lev[1:2]
datvarmelt2$subset3=factor(datvarmelt2$subset2,levels=lev1)



# create a theme for dot plots, which can be reused
theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


propclass=ggplot(datvarmelt2, aes(x = percent, y = subset3,color = variable, group = subset)) +
  geom_point() +
  theme_dotplot +
  scale_x_continuous(labels = scales::percent)+
    #scale_color_brewer(name="Model",palette="Dark2",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))
    scale_color_viridis(name="Model",discrete = TRUE, option = "C",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))


propclass = propclass +guides(size=FALSE)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 5.2, ymin = 4.8, fill = "grey", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 3.2, ymin = 1.8, fill = "palegreen", alpha = 0.2)+ ylab("gene sets")

#annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent),  ymax = "Dispensable", ymin = "Keystone", fill = "palegreen", alpha = 0.2)

ggsave("SplineULClasses.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineULClasses.pdf",propclass,width = 12, height = 7, units = "cm")



#### Same, removing 1 curve For Fig 3

datvarmelt2[datvarmelt2$variable!="X1.curve",] %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  mutate(freq = value / sumvar)  %>% 
  ungroup() -> datvarmelt2b

propclass=ggplot(datvarmelt2b, aes(x = freq, y = subset3,fill = variable)) +
  #geom_bar(position=position_dodge(),stat="identity")+
  geom_bar(stat="identity")+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))

ggsave("SplineULClassesBarplots.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineULClassesBarplots.pdf",propclass,width = 12, height = 7, units = "cm")



#### Same, divergence For Fig 5

datvarmelt2$divham="no"
datvarmelt2$divham[datvarmelt2$variable=="ham.2.teeth"|datvarmelt2$variable=="X2.teeth"]="ham"

datvarmelt2$divmus="no"
datvarmelt2$divmus[datvarmelt2$variable=="mus.2.teeth"|datvarmelt2$variable=="X2.teeth"]="mus"


datvarmelt2 %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divham=="ham")%>%
  mutate(type = "ham") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divmus=="mus")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2ham



datvarmelt2 %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divmus=="mus")%>%
  mutate(type = "mus") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divham=="ham")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2mus

datvarmelt2species=data.frame(rbind(datvarmelt2mus,datvarmelt2ham))

propclass=ggplot(datvarmelt2species, aes(x = freq, y = subset3,fill = type)) +
  geom_bar(position=position_dodge(),stat="identity")+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Ham 2curves","Mus 2curves"))

ggsave("SplineDivergentClassesBarplots.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineDivergentClassesBarplots.pdf",propclass,width = 12, height = 7, units = "cm")







### Proportion per signalling pathways

r=data.frame(lapply(names(pathwaylistspe),function(x){
  v=pathwaylistspe[[x]]
  n=getclassesmod3(list=v,suf=x,threshold=0.05)
  names(n)=c("subset","N","1 curve","mus 2 teeth","ham 2 teeth","2 teeth")
  n[is.na(n)]=0
  return(n)
}))
names(r)=names(pathwaylistspe)
r
library(forcats)
datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
row.names(datvar)=datvar$subset
datvarp=data.frame(t(r))
datvarp=datvarp[order(as.numeric(datvarp$N),decreasing=TRUE),]


datvarp=rbind(datvar[datvar$subset%in%c("Total","bite-it","Dispensable","Keystone","Pathways"),],datvarp)
datvarmeltp=melt(datvarp,id=c("subset","N"))
datvarmeltp$value=as.numeric(as.character(datvarmeltp$value))
datvarmelt2p <- datvarmeltp %>%
  group_by(subset) %>%
  mutate(percent = value / sum(value)) %>%
  ungroup()%>%
  mutate(subset2 = paste0(subset," (",N,")"))%>%
  mutate(subset2 = fct_reorder(subset2, as.numeric(as.character(N))))


lev= levels(datvarmelt2p$subset2)
lev1=lev
lev1[9]=lev[12]
lev1[10:12]=lev[9:11]

datvarmelt2p$subset3=factor(datvarmelt2p$subset2,levels=lev1)


propclass2=ggplot(datvarmelt2p, aes(x = percent, y =  subset3,color = variable, group = subset)) +
  geom_jitter(height = 0.1)+
  theme_dotplot +
  scale_x_continuous(labels = scales::percent)+
  scale_y_discrete(limits = rev(levels(subset)), labels=)+
  ylab("Gene subsets")+
  #scale_color_brewer(name="Model",palette="Dark2",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))
  scale_color_viridis(name="Model",discrete = TRUE, option = "C",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))

propclass2 = propclass2 +guides(size=FALSE)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 13.2, ymin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 11.2, ymin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 9.2, ymin = 0.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2p$percent)+0.1,  ymax = 9.2, ymin = 8.8, fill = "lightsteelblue", alpha = 0.2)
ggsave("SplineULClassesPathways.png",propclass2,width = 15, height = 10, units = "cm") 

ggsave("SplineULClassesPathways.pdf",propclass2,width = 15, height = 10, units = "cm") 
scp=propclass2




## barplots for UpLow diff


datvarmelt2p$divham="no"
datvarmelt2p$divham[datvarmelt2p$variable=="ham.2.teeth"|datvarmelt2p$variable=="X2.teeth"]="ham"

datvarmelt2p$divmus="no"
datvarmelt2p$divmus[datvarmelt2p$variable=="mus.2.teeth"|datvarmelt2p$variable=="X2.teeth"]="mus"


datvarmelt2p %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divham=="ham")%>%
  mutate(type = "ham") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divmus=="mus")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2pham

datvarmelt2p %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divmus=="mus")%>%
  mutate(type = "mus") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divham=="ham")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2pmus

datvarmelt2pspecies=data.frame(rbind(datvarmelt2pmus,datvarmelt2pham))

theme_barplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())



propclass=ggplot(datvarmelt2pspecies, aes(x = freq, y = subset3,fill = type)) + theme_barplot +
  geom_bar(width = 0.5,position=position_dodge(),stat="identity")+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Ham 2curves","Mus 2curves"))

ggsave("SplineDivergentClassesPatBarplots.png",propclass,width = 15, height = 7, units = "cm") 
ggsave("SplineDivergentClassesPatBarplots.pdf",propclass,width = 15, height = 7, units = "cm")



















### GO analysis 

lists=list(rownames(resmod3)[resmod3$mus=="Tooth"],rownames(resmod3)[resmod3$ham=="Tooth"],rownames(resmod3)[resmod3$model=="2 teeth"])
names(lists)=c("mus","ham","2teeth")

gostres <- gost(query = lists, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = rownames(resmod3), 
                numeric_ns = "", sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"), as_short_link = FALSE)


gostplot(gostres, capped = TRUE, interactive = TRUE)


gem <- gostres$result[,c("query","term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("model","GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem <- gem[,c("model","GO.ID", "Description", "FDR")]
write.table(gem, file = "gProfiler_gem_UpLoModels.txt", sep = "\t", quote = F, row.names = F)


# GO analysis  : Emaplot
gene.df <- bitr(rownames(resmod3), fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL","MGI"),
                OrgDb = org.Mm.eg.db)
df=merge(gene.df,resmod3,by.x="SYMBOL",by.y=0)


sapply(unique(df$model),function(cat){

ego <- enrichGO(gene          = df$ENSEMBL[df$model==cat],
                universe      = df$ENSEMBL,
                keyType       = 'ENSEMBL',
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
p=barplot(ego2, showCategory=50)
ggsave(paste0("Genes_",cat,"_barplot.pdf"),p,width=25,height=10,units="cm")

d <- godata('org.Mm.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
p=emapplot(ego2, layout="kk")
ggsave(paste0("Genes_",cat,"_emaplot.pdf"),p)

})





#### BCA ham sou 

load("ddsToothWhole.Rdata")

ddsToothWholeMus <- DESeqDataSetFromMatrix(countData = counts(ddsToothWhole)[,colData(ddsToothWhole)$espece=="mus"],
                              colData = colData(ddsToothWhole)[colData(ddsToothWhole)$espece=="mus",],
                              design = ~1)
colData(ddsToothWholeMus)$timerel=10*(as.numeric(colData(ddsToothWholeMus)$est_GAM)-14.5)/(18-14.5)
ddsToothWholeMus <- DESeq(ddsToothWholeMus)
pTwm=dudi.pca(t(counts(ddsToothWholeMus,norm=T)),scale=T,scann=F,n=10) 
bcasm=bca(pTwm,as.factor(as.character(colData(ddsToothWholeMus)$machoire)),nf=1,scannf=FALSE)
apply(pTwm$li,2, function(x){cor.test(x,bcasm$ls$CS1,meth="spearman")$estimate})
dfm=colData(ddsToothWholeMus)
dfm$bca=bcasm$ls$CS1
dfm1 = as.data.frame.matrix(dfm) %>% group_by(stade, rep) %>% filter(n() > 1) %>%
  mutate(Difference = bca - lag(bca))
bmus1=ggplot(dfm1,aes(x=as.numeric(stade),y=bca,shape=machoire,group=machoire))+geom_point(size=2)+xlab("stage")+ylab("BCA axis")+geom_smooth(se = FALSE)
bmus2=ggplot(dfm1[!is.na(dfm1$Difference),],aes(x=as.numeric(stade),y=Difference))+geom_point(size=2)+xlab("stage")+ylab("distance (Lower-Upper)")+geom_smooth(se = FALSE)





ddsToothWholeHam <- DESeqDataSetFromMatrix(countData = counts(ddsToothWhole)[,colData(ddsToothWhole)$espece=="ham"],
                                           colData = colData(ddsToothWhole)[colData(ddsToothWhole)$espece=="ham",],
                                           design = ~1)
colData(ddsToothWholeHam)$timerel=100*(as.numeric(colData(ddsToothWholeHam)$est_GAM)-12.2)/(14.5-12.2)
ddsToothWholeHam <- DESeq(ddsToothWholeHam)
pTh=dudi.pca(t(counts(ddsToothWholeHam,norm=T)),scale=T,scann=F,n=10) 
bcash=bca(pTh,as.factor(as.character(colData(ddsToothWholeHam)$machoire)),nf=1,scannf=FALSE)
dfh=colData(ddsToothWholeHam)
dfh$bca=-bcash$ls$CS1
dfh1 = as.data.frame.matrix(dfh) %>% group_by(stade, rep) %>% filter(n() > 1) %>%
  mutate(Difference = bca - lag(bca))
bham1=ggplot(dfh1,aes(x=as.numeric(stade),y=bca,shape=machoire,group=machoire))+geom_point(size=2)+xlab("stage")+ylab("BCA axis")+geom_smooth(se = FALSE)
bham2=ggplot(dfh1[!is.na(dfh1$Difference),],aes(x=as.numeric(stade),y=-Difference))+geom_point(size=2)+xlab("stage")+ylab("distance (Lower-Upper)")+geom_smooth(se = FALSE)



fig4 <- plot_grid(
  bmus2 + theme_light()+ theme(legend.position="none"),
  bham2 + theme_light()+ theme(legend.position="none"),
  align = 'vh',
  labels = c("A. Mus", "B. Ham"),
  hjust = -1,
  nrow = 1
)

ggsave("BCAstage2.png",fig4,width = 40, height = 20, units = "cm") 
ggsave("BCAstage2.pdf",fig4,width = 40, height = 20, units = "cm") 


