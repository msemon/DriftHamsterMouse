---
title: "BCA"
author: "marie"
date: "June 1st 2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r load, include=FALSE}
library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(readxl)
library(ade4)
library(dplyr)
library(forcats)
library(viridis)


## whole tooth data
metadataTot1=read.table(file="../Counts/metadataTot.txt",sep="\t")
CountTot=read.table(file="../Counts/CountTot.txt",sep="\t")
CountTot=CountTot[,metadataTot1$machoire=="mx"|metadataTot1$machoire=="md"]
metadataTot1=metadataTot1[metadataTot1$machoire=="mx"|metadataTot1$machoire=="md",]

metadataTot1$stade[metadataTot1$stade==12.2]=12.25
metadataTot1$ech=paste(metadataTot1$espece,metadataTot1$machoire,metadataTot1$stade,metadataTot1$rep,sep="_")
metadataTot1$id=row.names(metadataTot1)
metadataTot1=metadataTot1[,!names(metadataTot1)%in%c("stade")]



ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot,0),
                                        colData = metadataTot1,
                                        design = ~1)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="mus"]=10*(as.numeric(colData(ddsToothWhole)$stade[colData(ddsToothWhole)$espece=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="ham"]=10*(as.numeric(colData(ddsToothWhole)$stade[colData(ddsToothWhole)$espece=="ham"])-12)/(14.5-12)
ddsToothWhole <- DESeq(ddsToothWhole)




### lists of genes
biteit=read.table("liste_bite-it.csv",h=F)
keystone=read.csv("keystone_genes.csv",h=T,sep="\t")
dispensable=read.csv("dispensable_genes.csv",h=T,sep="\t")
pathway = read_excel("pathways-Margaux-corMS.xlsx",sheet=1)
pathway=pathway[,1:5]
pathwaylist=unique(pathway$Symbol)
pathwaylistspe=lapply(split(pathway,pathway$Pathway),function(x){unique(x$Symbol)})


ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot,0),
                                        colData = metadataTot1,
                                        design = ~ 1)
colData(ddsToothWhole)$espece=as.factor(colData(ddsToothWhole)$espece)     
colData(ddsToothWhole)$machoire=as.factor(colData(ddsToothWhole)$machoire)    
### HERE specify time to use
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="mus"]=10*(as.numeric(colData(ddsToothWhole)$est_GAM[colData(ddsToothWhole)$espece=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="ham"]=10*(as.numeric(colData(ddsToothWhole)$est_GAM[colData(ddsToothWhole)$espece=="ham"])-12.3)/(14.5-12.3)
ddsToothWhole <- DESeq(ddsToothWhole)
ddsToothWhole <- DESeq(ddsToothWhole)
### make DESeq object
save(file="ddsToothWhole.Rdata",ddsToothWhole)

```


```PCA MOUSE AND HAM SEPARATELY


```{r PCA tooth whole MUS}

ddsToothWholeMus <- DESeqDataSetFromMatrix(countData = round(CountTot[,metadataTot1$espece=="mus"&(metadataTot1$machoire=="md"|metadataTot1$machoire=="mx")],0),
                              colData = metadataTot1[metadataTot1$espece=="mus"&(metadataTot1$machoire=="md"|metadataTot1$machoire=="mx"),],
                              design = ~ 1)

                              
ddsToothWholeMus <- DESeq(ddsToothWholeMus)


ct=counts(ddsToothWholeMus,norm=T)
c=colData(ddsToothWholeMus)

samplelist=function(liste=biteit$V1,suf="bite-it")
{
lss=sapply(1:1000,function(s)
{
databite=ct[row.names(ct)%in%liste,]
s=sample(row.names(databite),replace=TRUE)
ps=dudi.pca(t(databite[s,]),scale=T,scann=F,n=10) 
psm=bca(ps,as.factor(as.character(c$machoire)),nf=1,scannf=FALSE)
maceffect=psm$eig[1]/sum(ps$eig)
return(c(maceffect*100))
})

databite=ct[row.names(ct)%in%liste,]
ps=dudi.pca(t(databite),scale=T,scann=F,n=10) 
psm=bca(ps,as.factor(as.character(c$machoire)),nf=1,scannf=FALSE)
maceffect=psm$eig[1]/sum(ps$eig)

lss=data.frame(t(lss))
q=quantile(lss,probs=c(0.025,0.975))[1:2]
return(c(maceffect,q[1],q[2],suf))
}


s1=samplelist(liste=biteit$V1,suf="bite-it")
s2=samplelist(liste=dispensable$Gene..name,suf="Dispensable")
s3=samplelist(liste=keystone$Gene.name,suf="Keystone")
s6=samplelist(liste=row.names(CountTot),suf="Total")
s7=samplelist(liste=pathwaylist,suf="Pathways")



datvar=data.frame(rbind(s1,s2,s3,s6,s7))
names(datvar)=c("value","lower","upper","subset")
datvar$value=as.numeric(as.character(datvar$value))
datvar$lower=as.numeric(as.character(datvar$lower))/100
datvar$upper=as.numeric(as.character(datvar$upper))/100
datvar$subset=factor(datvar$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))
datvarmus=datvar


# create a theme for dot plots, which can be reused
theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


propclass=ggplot(datvar, aes(x = value, y = subset, group = subset)) +
  geom_point() +
  theme_dotplot +
  scale_x_continuous(labels = scales::percent)+
   geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                 position=position_dodge(.9)) 
   # scale_color_viridis(name="Model",discrete = TRUE, option = "C",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))

ggsave("variancetooth_Mus.png",propclass,width = 20, height = 15, units = "cm") 


````


```{r PCA tooth whole HAM}

ddsToothWholeHam <- DESeqDataSetFromMatrix(countData = round(CountTot[,metadataTot1$espece=="ham"&(metadataTot1$machoire=="md"|metadataTot1$machoire=="mx")],0),
                              colData = metadataTot1[metadataTot1$espece=="ham"&(metadataTot1$machoire=="md"|metadataTot1$machoire=="mx"),],
                              design = ~ 1)
ddsToothWholeHam <- DESeq(ddsToothWholeHam)


ct=counts(ddsToothWholeHam,norm=T)
c=colData(ddsToothWholeHam)

samplelist=function(liste=biteit$V1,suf="bite-it")
{
lss=sapply(1:1000,function(s)
{
databite=ct[row.names(ct)%in%liste,]
s=sample(row.names(databite),replace=TRUE)
ps=dudi.pca(t(databite[s,]),scale=T,scann=F,n=10) 
psm=bca(ps,as.factor(as.character(c$machoire)),nf=1,scannf=FALSE)
maceffect=psm$eig[1]/sum(ps$eig)
return(c(maceffect*100))
})

databite=ct[row.names(ct)%in%liste,]
ps=dudi.pca(t(databite),scale=T,scann=F,n=10) 
psm=bca(ps,as.factor(as.character(c$machoire)),nf=1,scannf=FALSE)
maceffect=psm$eig[1]/sum(ps$eig)

lss=data.frame(t(lss))
q=quantile(lss,probs=c(0.025,0.975))[1:2]
return(c(maceffect,q[1],q[2],suf))
}


s1=samplelist(liste=biteit$V1,suf="bite-it")
s2=samplelist(liste=dispensable$Gene..name,suf="Dispensable")
s3=samplelist(liste=keystone$Gene.name,suf="Keystone")
s6=samplelist(liste=row.names(CountTot),suf="Total")
s7=samplelist(liste=pathwaylist,suf="Pathways")



datvar=data.frame(rbind(s1,s2,s3,s6,s7))
names(datvar)=c("value","lower","upper","subset")
datvar$value=as.numeric(as.character(datvar$value))
datvar$lower=as.numeric(as.character(datvar$lower))/100
datvar$upper=as.numeric(as.character(datvar$upper))/100
datvar$subset=factor(datvar$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))
datvarham=datvar


# create a theme for dot plots, which can be reused
theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


propclass=ggplot(datvar, aes(x = value, y = subset, group = subset)) +
  geom_point() +
  theme_dotplot +
  scale_x_continuous(labels = scales::percent)+
   geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                 position=position_dodge(.9)) 
   # scale_color_viridis(name="Model",discrete = TRUE, option = "C",labels = c("1 curve", "Mus 2 curves", "Ham 2 curves","MusHam 2 curves"))

ggsave("variancetooth_Ham.png",propclass,width = 20, height = 15, units = "cm") 


````


```{r Fig Ham Mus}


datvarham$species="ham"
datvarmus$species="mus"
datvar=data.frame(rbind(datvarham,datvarmus))

propclass=ggplot(datvar, aes(x = value, y = subset, fill=species)) +
  theme_dotplot +
  geom_bar(width = 0.5,position=position_dodge(),stat="identity")+
  scale_x_continuous(labels = scales::percent)+
   geom_errorbar(width = 0.5,aes(xmin=lower, xmax=upper),
                 position=position_dodge()) +
   scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Hamster", "Mouse"))

ggsave("BCA_ham_mus.png",propclass,width = 20, height = 15, units = "cm") 
ggsave("BCA_ham_mus.pdf",propclass,width = 20, height = 15, units = "cm") 

propclass=ggplot(datvar, aes(x = value, y = subset, color=species)) +
  theme_dotplot +
  geom_point(position=position_dodge(width = 0.5),stat="identity",size=2)+
  scale_x_continuous(labels = scales::percent)+
   geom_errorbar(width = 0.5,aes(xmin=lower, xmax=upper),
                 position=position_dodge()) +
   scale_color_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Hamster", "Mouse"))

ggsave("BCA_ham_mus_colo.png",propclass,width = 20, height = 15, units = "cm") 
ggsave("BCA_ham_mus_colo.pdf",propclass,width = 20, height = 15, units = "cm") 



````
