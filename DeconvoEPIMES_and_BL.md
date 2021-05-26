---
title: "Deconvolutions"
author: "marie"
date: "May 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Data Preparation

```{r prep}

library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(ade4)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(DeconRNASeq)
```


```{r prep}

metadataTot1=read.table(file="metadataTot_FINAL.txt")
CountTot=read.table(file="CountTot.txt")

# to remove kidney samples used as controls
ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot[,-c(grep("kidney",metadataTot1$machoire))],0),
                              colData = metadataTot1[-c(grep("kidney",metadataTot1$machoire)),],
                              design = ~ 1)

colData(ddsToothWhole)$stade=as.numeric(as.character(colData(ddsToothWhole)$stade))
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="mus"]=10*(as.numeric(colData(ddsToothWhole)$stade[colData(ddsToothWhole)$espece=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$espece=="ham"]=10*(as.numeric(colData(ddsToothWhole)$stade[colData(ddsToothWhole)$espece=="ham"])-12)/(14.5-12)
colData(ddsToothWhole)$machoire=as.factor(colData(ddsToothWhole)$machoire)                              
                              
ddsToothWhole <- DESeq(ddsToothWhole)
vstT=data.frame(assay(vst(ddsToothWhole)))
bm=data.frame(counts(ddsToothWhole,norm=T))



## time estimates from GAM model
mdhamtime=read.table("predictions_rnaseq_cuspmdH_all.txt",h=T)
mxhamtime=read.table("predictions_rnaseq_cuspmxH_all.txt",h=T)
mdmustime=read.table("predictions_rnaseq_cuspmdM_all.txt",h=T)
mxmustime=read.table("predictions_rnaseq_cuspmxM_all.txt",h=T)

times=data.frame(rbind(mdhamtime,mxhamtime,mdmustime,mxmustime))
times$echantillon[grep("mus",times$echantillon)]=paste0(times$echantillon[grep("mus",times$echantillon)],"W")

metadataTot2=merge(metadataTot1,times,by.x=0,by.y="echantillon",all.x=T,no.dups=TRUE)
row.names(metadataTot2)=metadataTot2$Row.names
metadataTot2=metadataTot2[row.names(metadataTot1),]

```


# Epithelium/mesenchyme specific genes

```{r DEepimes}

ddsToothPure <- DESeqDataSetFromMatrix(countData = round(CountTot[,c(grep("mes",metadataTot2$machoire),grep("epi",metadataTot2$machoire))],0),
                              colData = metadataTot2[c(grep("mes",metadataTot2$machoire),grep("epi",metadataTot2$machoire)),],
                              design = ~ 1)


colData(ddsToothPure)$tissue="mes"
colData(ddsToothPure)$tissue[grep("epi",colData(ddsToothPure)$machoire)]="epi"

ddsToothPure <- DESeqDataSetFromMatrix(countData = round(CountTot[,c(grep("mes",metadataTot2$machoire),grep("epi",metadataTot2$machoire))],0),
                              colData = colData(ddsToothPure),
                              design = ~ tissue)
ddsToothPure <- DESeq(ddsToothPure)                              
respt=results(ddsToothPure)

table(respt$padj<0.05) 

metadataTot3=metadataTot1[c(grep("mes",metadataTot2$machoire),grep("epi",metadataTot2$machoire)),]
metadataTot3$tis[grep("mes",metadataTot3$machoire)]="mes"
metadataTot3$tis[grep("epi",metadataTot3$machoire)]="epi"

basemean_counts=as.data.frame(counts(ddsToothPure,norm=T))



```


```{r deconvolutions epithelium/mesenchym}

fdeconvo=function(seuilL2FC=5,seuilPV=0.05,nboot=100)
{
signaturesnames=respt[!is.na(respt$padj)&abs(respt$log2FoldChange)>seuilL2FC&respt$padj<seuilPV,]
e=as.vector(apply(basemean_counts[row.names(basemean_counts)%in%row.names(signaturesnames),metadataTot3$tis=="epi"],1,mean))
m=as.vector(apply(basemean_counts[row.names(basemean_counts)%in%row.names(signaturesnames),metadataTot3$tis=="mes"],1,mean))
signaturesBM=as.data.frame(cbind(epi=e,mes=m))
row.names(signaturesBM)=row.names(signaturesnames)


# control
basemean_counts_tot=as.data.frame(counts(ddsToothPure,norm=T))
d=DeconRNASeq(basemean_counts_tot, signaturesBM, checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)
deconv_verif=as.data.frame(cbind(colData(ddsToothPure),d$out.all))
deconv_verif$tooth="lower"
deconv_verif$tooth[grep("mx",deconv_verif$machoire)]="upper"
deconv_verif$tissue="mes"
deconv_verif$tissue[grep("epi",deconv_verif$machoire)]="epi"

p=ggplot(deconv_verif, aes(x=tooth, y=mes, group=stade.x)) + 
  facet_grid(espece ~ tissue)+theme_bw()+
   geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Tooth type",y="% mesenchyme (deconvolutions)")
ggsave(filename=paste("deconv_control_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),p)


basemean_counts_tot=as.data.frame(counts(ddsToothWhole,norm=T))
metadata_counts_tot=merge(colData(ddsToothWhole),times,by.x=0,by.y="echantillon",all.x=T,no.dups=TRUE)
row.names(metadata_counts_tot)=metadata_counts_tot$Row.names
metadata_counts_tot=metadata_counts_tot[row.names(colData(ddsToothWhole)),]


#bootstraps
bootDec=function(prop=0.5)
{
  s=sample(1:nrow(signaturesBM),round(nrow(signaturesBM)*prop,0))
  d=DeconRNASeq(basemean_counts_tot, signaturesBM[s,], checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)$out.all
  return(d)
}

z=lapply(1:nboot,function(x){bootDec(prop=0.5)})
zz=data.frame(lapply(z,function(x){x[,"mes"]}))
zzz=apply(zz,1,function(s){quantile(s,p=c(0,0.25,0.5,0.75,1))})

deconv_tab=as.data.frame(cbind(metadata_counts_tot,mini=zzz[2,],median=zzz[3,],maxi=zzz[4,]))
deconv_tab$tooth="lower"
deconv_tab$tooth[grep("mx",deconv_tab$machoire)]="upper"

deconv_tab$type="part"
deconv_tab$type[deconv_tab$machoire=="md"]="whole"
deconv_tab$type[deconv_tab$machoire=="mx"]="whole"


deconv_tab$timerelGAM[deconv_tab$espece=="mus"]=10*(deconv_tab$est_GAM[deconv_tab$espece=="mus"]-14.5)/(18-14.5)
deconv_tab$timerelGAM[deconv_tab$espece=="ham"]=10*(deconv_tab$est_GAM[deconv_tab$espece=="ham"]-12)/(14.5-12)



p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=est_GAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=est_GAM, xend=est_GAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% mesenchyme")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_mes_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)
                               
 
 

p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerelGAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerelGAM, xend=timerelGAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% mesenchyme")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_mes_deconRNAseq_timeGAM_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)
                               
 
 
 
 
p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerel, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerel, xend=timerel,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% mesenchyme")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))                              
                               
save_plot(filename=paste("deconv_total_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)

write.table(file=paste0("deconv_",seuilL2FC,"pv",seuilPV,".txt"),deconv_tab,sep="\t",quote=F)

return(list(p,deconv_tab))
}



boot5=fdeconvo(seuilL2FC=5,seuilPV=0.05,nboot=1000)
boot7=fdeconvo(seuilL2FC=7,seuilPV=0.05,nboot=1000)
boot3=fdeconvo(seuilL2FC=3,seuilPV=0.05,nboot=1000)


deconv_tab=boot3[[2]]
p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerelGAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerelGAM, xend=timerelGAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% mesenchyme")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_mes_deconRNAseq_timeGAM_logFC",seuilL2FC,"pv",seuilPV,".pdf",sep=""),ncol=2,nrow=1,p)
                               

```




Deconvolution bucco-lingual


# Bucco/lingual specific genes

```{r bl markers}

ddsToothPure <- DESeqDataSetFromMatrix(countData = round(CountTot[,c(grep("_B",metadataTot2$machoire),grep("_L",metadataTot2$machoire))],0),
                              colData = metadataTot2[c(grep("_B",metadataTot2$machoire),grep("_L",metadataTot2$machoire)),],
                              design = ~ 1)


colData(ddsToothPure)$tissue="_B"
colData(ddsToothPure)$tissue[grep("_L",colData(ddsToothPure)$machoire)]="_L"

ddsToothPure <- DESeqDataSetFromMatrix(countData = round(CountTot[,c(grep("_B",metadataTot2$machoire),grep("_L",metadataTot2$machoire))],0),
                              colData = colData(ddsToothPure),
                              design = ~ tissue)
ddsToothPure <- DESeq(ddsToothPure)                              
respt=results(ddsToothPure)

table(respt$padj<0.05) 

metadataTot3=metadataTot1[c(grep("_B",metadataTot2$machoire),grep("_L",metadataTot2$machoire)),]
metadataTot3$tis[grep("_B",metadataTot3$machoire)]="B"
metadataTot3$tis[grep("_L",metadataTot3$machoire)]="L"

basemean_counts=as.data.frame(counts(ddsToothPure,norm=T))
```


```{r bucco/lingual deconvolutions}

fdeconvoBL=function(seuilL2FC=5,seuilPV=0.05,nboot=100)
{
signaturesnames=respt[!is.na(respt$padj)&abs(respt$log2FoldChange)>seuilL2FC&respt$padj<seuilPV,]
b=as.vector(apply(basemean_counts[row.names(basemean_counts)%in%row.names(signaturesnames),metadataTot3$tis=="B"],1,mean))
l=as.vector(apply(basemean_counts[row.names(basemean_counts)%in%row.names(signaturesnames),metadataTot3$tis=="L"],1,mean))
signaturesBM=as.data.frame(cbind(buc=b,ling=l))
row.names(signaturesBM)=row.names(signaturesnames)


# control
basemean_counts_tot=as.data.frame(counts(ddsToothPure,norm=T))
d=DeconRNASeq(basemean_counts_tot, signaturesBM, checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)
deconv_verif=as.data.frame(cbind(colData(ddsToothPure),d$out.all))
deconv_verif$tooth="lower"
deconv_verif$tooth[grep("mx",deconv_verif$machoire)]="upper"
deconv_verif$tissue="B"
deconv_verif$tissue[grep("L",deconv_verif$machoire)]="L"

p=ggplot(deconv_verif, aes(x=tooth, y=buc, group=stade.x)) + 
  facet_grid(espece ~ tissue)+theme_bw()+
   geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Tooth type",y="% buccal (deconvolutions)")
ggsave(filename=paste("deconvBL_control_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),p)


basemean_counts_tot=as.data.frame(counts(ddsToothWhole,norm=T))
metadata_counts_tot=merge(colData(ddsToothWhole),times,by.x=0,by.y="echantillon",all.x=T,no.dups=TRUE)
row.names(metadata_counts_tot)=metadata_counts_tot$Row.names
metadata_counts_tot=metadata_counts_tot[row.names(colData(ddsToothWhole)),]


#bootstraps
bootDec=function(prop=0.5)
{
  s=sample(1:nrow(signaturesBM),round(nrow(signaturesBM)*prop,0))
  d=DeconRNASeq(basemean_counts_tot, signaturesBM[s,], checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)$out.all
  return(d)
}

z=lapply(1:nboot,function(x){bootDec(prop=0.5)})
zz=data.frame(lapply(z,function(x){x[,"buc"]}))
zzz=apply(zz,1,function(s){quantile(s,p=c(0,0.25,0.5,0.75,1))})

deconv_tab=as.data.frame(cbind(metadata_counts_tot,mini=zzz[2,],median=zzz[3,],maxi=zzz[4,]))
deconv_tab$tooth="lower"
deconv_tab$tooth[grep("mx",deconv_tab$machoire)]="upper"

deconv_tab$type="part"
deconv_tab$type[deconv_tab$machoire=="md"]="whole"
deconv_tab$type[deconv_tab$machoire=="mx"]="whole"


deconv_tab$timerelGAM[deconv_tab$espece=="mus"]=10*(deconv_tab$est_GAM[deconv_tab$espece=="mus"]-14.5)/(18-14.5)
deconv_tab$timerelGAM[deconv_tab$espece=="ham"]=10*(deconv_tab$est_GAM[deconv_tab$espece=="ham"]-12)/(14.5-12)



p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=est_GAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=est_GAM, xend=est_GAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% buccal")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_buc_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)
                               
 
 

p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerelGAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerelGAM, xend=timerelGAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% buccal")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_buc_deconRNAseq_timeGAM_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)
                               
 
 
 
 
p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerel, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerel, xend=timerel,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% buccal")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))                              
                               
save_plot(filename=paste("deconv_buc_total_deconRNAseq_logFC",seuilL2FC,"pv",seuilPV,".png",sep=""),ncol=2,nrow=1,p)

write.table(file=paste0("deconv_buc_",seuilL2FC,"pv",seuilPV,".txt"),deconv_tab,sep="\t",quote=F)

return(list(p,deconv_tab))
}

bootBL1=fdeconvoBL(seuilL2FC=1,seuilPV=0.05,nboot=1000)
bootBL0.5=fdeconvoBL(seuilL2FC=0.5,seuilPV=0.05,nboot=1000)
bootBL0=fdeconvoBL(seuilL2FC=0,seuilPV=0.05,nboot=1000)


seuilL2FC=1
deconv_tab=bootBL1[[2]]
p=ggplot(deconv_tab[deconv_tab$type=="whole",], aes(x=timerelGAM, y=median, color=tooth)) + 
  facet_grid(. ~ espece, scales = "free_x")+theme_bw()+
   geom_point(size=2,alpha = 0.8) + theme_bw(base_size = 22)+
   ylim(0,1)+
   geom_segment(data=deconv_tab[deconv_tab$type=="whole",],aes(x=timerelGAM, xend=timerelGAM,y=mini,yend=maxi)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
         labs(x="Dev. time",y="% mesenchyme")+
         scale_colour_manual(name="Tooth",
                      values=c("upper"="black",
                               "lower"="grey"))
save_plot(filename=paste("deconv_bl_deconRNAseq_timeGAM_logFC",seuilL2FC,"pv",seuilPV,".pdf",sep=""),ncol=2,nrow=1,p)
                               

```


