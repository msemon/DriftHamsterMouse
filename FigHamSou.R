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


resmod2=data.frame(mod2.md.padj=res_md_vs_sp[[2]]$padj,
                   mod2.md.dist=res_md_vs_sp[[2]]$dist,
                   mod2.mx.padj=res_mx_vs_sp[[2]]$padj,
                   mod2.mx.dist=res_mx_vs_sp[[2]]$dist)
row.names(resmod2)=row.names(res_md_vs_sp[[2]])

resmod2$md="1 curve"
resmod2$md[resmod2$mod2.md.padj<0.05]="Species"
resmod2$mx="1 curve"
resmod2$mx[resmod2$mod2.mx.padj<0.05]="Species"
resmod2$model="1 curve"
resmod2$model[resmod2$md=="Species"&resmod2$mx=="Species"]="LowerUpper"
resmod2$model[resmod2$md=="Species"&resmod2$mx=="1 curve"]="Lower"
resmod2$model[resmod2$mx=="Species"&resmod2$md=="1 curve"]="Upper"


### 

pathway = read_excel("../pathways-Margaux-corMS.xlsx",sheet=1)
pathway=pathway[,1:5]

pathwaylist=unique(pathway$Symbol)
pathwaylistspe=lapply(split(pathway,pathway$Pathway),function(x){unique(x$Symbol)})


 
 getclassesmod2=function(list=biteit$V1,suf="bite-it",threshold=0.05)
 {
   subset=resmod2[row.names(resmod2)%in%list,]
   out=c(table(subset$model),nrow(subset),suf)
   names(out)=c(names(table(subset$model)),"N","subset")
   out=out[c("subset","N","1 curve","Lower","Upper","LowerUpper")]
   return(out)
 }
 

n1=getclassesmod2(list=biteit$V1,suf="bite-it")
n2=getclassesmod2(list=dispensable$Gene..name,suf="Dispensable")
n3=getclassesmod2(list=keystone$Gene.name,suf="Keystone")
n6=getclassesmod2(list=row.names(resmod3),suf="Total")
n7=getclassesmod2(list=pathwaylist,suf="Pathways")


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
  scale_color_viridis(name="Model",discrete = TRUE, option = "D",labels = c("1 curve","Lower","Upper","LowerUpper"))

    #scale_color_brewer(name="Model",palette="Dark2",labels = c("1 curve","Lower","Upper","LowerUpper"))

propclass = propclass +guides(size=FALSE)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 5.2, ymin = 4.8, fill = "grey", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 3.2, ymin = 1.8, fill = "palegreen", alpha = 0.2)+ ylab("gene sets")

#annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent),  ymax = "Dispensable", ymin = "Keystone", fill = "palegreen", alpha = 0.2)

ggsave("SplineHamMusClasses.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineHamMusClasses.pdf",propclass,width = 12, height = 7, units = "cm")


### Proportion per signalling pathways

r=data.frame(lapply(names(pathwaylistspe),function(x){
  v=pathwaylistspe[[x]]
  n=getclassesmod2(list=v,suf=x,threshold=0.05)
  names(n)=c("subset","N","1 curve","Lower","Upper","LowerUpper")
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
  scale_color_viridis(name="Model",discrete = TRUE, option = "D",labels = c("1 curve","Lower","Upper","LowerUpper"))
  #scale_color_brewer(name="Model",palette="Dark2",labels = c("1 curve","Lower","Upper","LowerUpper"))

propclass2 = propclass2 +guides(size=FALSE)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 13.2, ymin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 11.2, ymin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 9.2, ymin = 0.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2p$percent)+0.1,  ymax = 9.2, ymin = 8.8, fill = "lightsteelblue", alpha = 0.2)
ggsave("SplineMusHamClassesPathways.png",propclass2,width = 15, height = 10, units = "cm") 

ggsave("SplineMusHamClassesPathways.pdf",propclass2,width = 15, height = 10, units = "cm") 
scp=propclass2


########
## barplots for UpLow diff


datvarmelt2p$divlo="no"
datvarmelt2p$divlo[datvarmelt2p$variable=="Lower"|datvarmelt2p$variable=="LowerUpper"]="Lower"

datvarmelt2p$divup="no"
datvarmelt2p$divup[datvarmelt2p$variable=="Upper"|datvarmelt2p$variable=="LowerUpper"]="Upper"


datvarmelt2p %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divlo=="Lower")%>%
  mutate(type = "Lower") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divup=="Upper")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2pLow

datvarmelt2p %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  subset(divup=="Upper")%>%
  mutate(type = "Upper") %>%
  mutate(sumdiv = sum(value)) %>%
  subset(divlo=="Lower")%>%
  mutate(freq = sumdiv / sumvar) %>%
  ungroup() -> datvarmelt2pUp

datvarmelt2ptooth=data.frame(rbind(datvarmelt2pLow,datvarmelt2pUp))

theme_barplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


propclass=ggplot(datvarmelt2ptooth, aes(x = freq, y = subset3,fill = type)) + theme_barplot +
  geom_bar(width = 0.5,position=position_dodge(),stat="identity")+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "D",labels = c("Lower 2curves","Upper 2curves"))

ggsave("SplineMusHamClassesPatBarplots.png",propclass,width = 15, height = 7, units = "cm") 
ggsave("SplineMusHamClassesPatBarplots.pdf",propclass,width = 15, height = 7, units = "cm")




#### proportion of up&lo common for Fig 5

datvarmelt2p[datvarmelt2p$variable!="X1.curve",] %>%  
  group_by(subset) %>%
  mutate(sumvar = sum(value)) %>%
  mutate(freq = value / sumvar)  %>% 
  ungroup() -> datvarmelt2b

propclass=ggplot(datvarmelt2b, aes(x = freq, y = subset3,fill = variable)) +
  geom_bar(stat="identity")+ theme_barplot +
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "D")

ggsave("SplineClassesBarplots.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineClassesBarplots.pdf",propclass,width = 12, height = 7, units = "cm")












### GO analysis 


lists=list(rownames(resmod2)[resmod2$md=="Species"],rownames(resmod2)[resmod2$mx=="Species"],rownames(resmod2)[resmod2$model=="LowerUpper"])
names(lists)=c("Lower","Upper","LowerUpper")

gostres <- gost(query = lists, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"), as_short_link = FALSE)


gostplot(gostres, capped = TRUE, interactive = TRUE)


gem <- gostres$result[,c("query","term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("model","GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem <- gem[,c("model","GO.ID", "Description", "FDR")]
write.table(gem, file = "gProfiler_gem_MusHamModels.txt", sep = "\t", quote = F, row.names = F)



### GO analysis 

lists=list(rownames(resmod2)[resmod2$md=="Species"],rownames(resmod2)[resmod2$mx=="Species"],rownames(resmod2)[resmod2$model=="LowerUpper"])
names(lists)=c("md","mx","2teeth")

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
write.table(gem, file = "gProfiler_gem_MusHamModels.txt", sep = "\t", quote = F, row.names = F)


# GO analysis  : Emaplot
gene.df <- bitr(rownames(resmod2), fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL","MGI"),
                OrgDb = org.Mm.eg.db)
df=merge(gene.df,resmod2,by.x="SYMBOL",by.y=0)


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

