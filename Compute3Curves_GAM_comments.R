##### plot full models for a selection of genes and COMPUTE models with 3 curves (the top part of the code is redundant with Compute4Curves_GAM but adds the 3 curves models and some plots, plus enrichments

################################################################### load libraries
library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(readxl)
library(forcats)
library(dplyr)
library(viridis)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(readxl)
library(forcats)
library(dplyr)
library(viridis)
library(sjPlot)
library(org.Mm.eg.db)
library(gprofiler2)
library(clusterProfiler)
library(GOSemSim)
library(viridis)
library(readxl)
library(enrichplot)
library(simplifyEnrichment)

################################################################### load data

# metadata
metadataTot1=read.table(file="data/metadataTot.txt",sep="\t")
names(metadataTot1)=c("jaw","stage","species","rep","file")
# raw counts for 64 samples
CountTot=read.table(file="data/CountTot.txt",sep="\t")
CountTot=CountTot[,metadataTot1$jaw=="mx"|metadataTot1$jaw=="md"]
metadataTot1=metadataTot1[metadataTot1$jaw=="mx"|metadataTot1$jaw=="md",]

metadataTot1$stage[metadataTot1$stage==12.2]=12.25
metadataTot1$ech=paste(metadataTot1$species,metadataTot1$jaw,metadataTot1$stage,metadataTot1$rep,sep="_")
metadataTot1$id=row.names(metadataTot1)
metadataTot1=metadataTot1[,!names(metadataTot1)%in%c("stage")]


## time estimates from GAM model
mdhamtime=read.table("data/predictions_rnaseq_cuspmdH_all.txt",h=T)
mxhamtime=read.table("data/predictions_rnaseq_cuspmxH_all.txt",h=T)
mdmustime=read.table("data/predictions_rnaseq_cuspmdM_all.txt",h=T)
mxmustime=read.table("data/predictions_rnaseq_cuspmxM_all.txt",h=T)
times=data.frame(rbind(mdhamtime,mxhamtime,mdmustime,mxmustime))
names(times)=c("samples","stage","replicate","weight","fit_boxCox","lwr_boxCox","upr_boxCox","fit_log","lwr_log","upr_log","est_GAM","lower_GAM","upper_GAM")
times$sample[grep("mus",times$sample)]=paste0(times$sample[grep("mus",times$sample)],"W")

# add these estimated development times to colData and reorder samples in countdata
metadataTot2=merge(metadataTot1,times,by.x=0,by.y="sample",all.x=T,no.dups=TRUE)
CountTot2=CountTot[,metadataTot2$id]
row.names(metadataTot2)=metadataTot2$id

ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot2,0),
                                        colData = metadataTot2,
                                        design = ~1)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="mus"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="ham"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="ham"])-12)/(14.5-12)
ddsToothWhole <- DESeq(ddsToothWhole)

### lists of genes of interest 
biteit=read.table("data/liste_bite-it.csv",h=F)
keystone=read.csv("data/keystone_genes.csv",h=T,sep="\t")
dispensable=read.csv("data/dispensable_genes.csv",h=T,sep="\t")
pathway = read_excel("data/pathways-Margaux-corMS.xlsx",sheet=1)
pathway=pathway[,1:5]
pathwaylist=unique(pathway$Symbol)
pathwaylistspe=lapply(split(pathway,pathway$Pathway),function(x){unique(x$Symbol)})


ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot2,0),
                                        colData = metadataTot2,
                                        design = ~ 1)
colData(ddsToothWhole)$species=as.factor(colData(ddsToothWhole)$species)     
colData(ddsToothWhole)$jaw=as.factor(colData(ddsToothWhole)$jaw) 

### HERE specify time to use (could be read days, or PCA etc). We choose GAM here.
# make boundaries comparable between species (cap/bell)

colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="mus"]=10*(as.numeric(colData(ddsToothWhole)$est_GAM[colData(ddsToothWhole)$species=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="ham"]=10*(as.numeric(colData(ddsToothWhole)$est_GAM[colData(ddsToothWhole)$species=="ham"])-12.3)/(14.5-12.3)
ddsToothWhole <- DESeq(ddsToothWhole)


################################################################### Fit the models
# 2 boundary knots, polynomials degree 3
# DEsplines Time:(species,jaw)}

spline_design_matrix <- function(time, intercept = F) {
  tmin <- min(time)
  tmax <- max(time)
  X <- bs(time,
          Boundary.knots = c(tmin,tmax),
          #knots=c(0,3,6)
          intercept = intercept,
          degree = 3)
  colnames(X) <- paste("bs_", colnames(X), sep = "")
  X
}


# to get comparable stages in both species we remove 2 younger samples in hamster
df1=colData(ddsToothWhole)[colData(ddsToothWhole)$stage>=12,]
SDM <- spline_design_matrix(df1$timerel, intercept = F)
spline_factors <- colnames(SDM)

df1=cbind(df1,SDM)      
df1$jaw=factor(df1$jaw)
# factor coding upper Mouse versus the 3 other teeth (for model with tooth specific profiles in mouse)
df1$isupMus="no"
df1$isupMus[df1$species=="mus"&df1$jaw=="mx"]="yes"

# factor coding upper Mouse versus the 3 other teeth (for model with tooth specific profiles in hamster)
df1$isupHam="no"
df1$isupHam[df1$species=="ham"&df1$jaw=="mx"]="yes"

# complete model with intercept per species and interaction species/jaw
design.species.jaw <- as.formula(paste(" ~species  + (",
                                       paste(spline_factors,collapse=" + "),
                                       "):paste(species,jaw) + 0 "))

dds <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.species.jaw  ) 
dds <- DESeq(dds)
res.diff.lev_species <- results(dds,contrast=c("species","ham","mus"))

# species model with intercept per species and one curve per species
design.species <- as.formula(paste("~ species + (",
                                   paste(spline_factors,collapse=" + "),
                                   "):species + 0 "))
# test if full model is bettern than species model
dds.diff.sp <- DESeq(dds, test="LRT", reduced=design.species)
res.diff.pat_species <- results(dds.diff.sp)


# jaw model with intercept per species and one curve per jaw
design.jaw <- as.formula(paste("~ species + (",
                               paste(spline_factors,collapse=" + "),
                               "):jaw + 0 "))
# test if full model is bettern than jaw model
dds.diff.jaw <- DESeq(dds, test="LRT", reduced=design.jaw)
res.diff.pat_jaw <- results(dds.diff.jaw)


# simple model with intercept per species and same curve for all teeth
design.simple <- as.formula(paste("~  species + (",
                                  paste(spline_factors,collapse=" + "),
                                  ") + 0 "))
# test if species is better than simple model
dds.diff.simple <- DESeq(dds, test="LRT", reduced=design.simple)
res.diff.pat_simple <- results(dds.diff.simple)


# model with up/low difference in mouse, but not in hamster
design.3teethMusUp <- as.formula(paste("~ species + (",
                                       paste(spline_factors,collapse=" + "),
                                       "):paste(species,isupMus) + 0 "))
dds.3teethMusUp <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethMusUp  ) 
dds.3teethMusUp <- DESeq(dds.3teethMusUp, test="LRT", reduced=design.species)
res.diff.3teethMusUp <- results(dds.3teethMusUp)


# model with up/low difference in hamster, but not in mouse
design.3teethHamUp <- as.formula(paste("~ species + (",
                                       paste(spline_factors,collapse=" + "),
                                       "):paste(species,isupHam) + 0 "))
dds.3teethHamUp <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethHamUp  ) 
dds.3teethHamUp <- DESeq(dds.3teethHamUp, test="LRT", reduced=design.species)
res.diff.3teethHamUp <- results(dds.3teethHamUp)



## coeffs of the complete model (stored in dds object)

spline_coefs <- function(sp) {
  sapply(spline_factors,
         function(f) results(dds, name=paste0(f,".paste.species..jaw.",sp))[,2])
}
sc.ham.md <- spline_coefs("ham.md")
sc.ham.mx <- spline_coefs("ham.mx")
sc.mus.md <- spline_coefs("mus.md")
sc.mus.mx <- spline_coefs("mus.mx")
row.names(sc.ham.md)=row.names(dds)
row.names(sc.ham.mx)=row.names(dds)
row.names(sc.mus.md)=row.names(dds)
row.names(sc.mus.mx)=row.names(dds)


sc.mus <- results(dds, name=paste0("species","mus"))[,2]
sc.ham <- results(dds, name=paste0("species","ham"))[,2]
names(sc.mus)=row.names(dds)
names(sc.ham)=row.names(dds)

## coeffs of the species model

dds.sp <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.species  ) 
dds.sp <- DESeq(dds.sp)

spline_coefs.species <- function(sp) {
  sapply(spline_factors,
         function(f) results(dds.sp, name=paste0("species",sp,".",f))[,2])
}

sc.sp.mus <- spline_coefs.species("mus")
sc.sp.ham <- spline_coefs.species("ham")
row.names(sc.sp.mus)=row.names(dds.sp)
row.names(sc.sp.ham)=row.names(dds.sp)

sc.sp.fmus <- results(dds.sp, name=paste0("species","mus"))[,2]
sc.sp.fham <- results(dds.sp, name=paste0("species","ham"))[,2]
names(sc.sp.fmus)=row.names(dds)
names(sc.sp.fham)=row.names(dds)


## coeffs of the jaw model

dds.jaw <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.jaw  ) 
dds.jaw <- DESeq(dds.jaw)

spline_coefs <- function(sp) {
  sapply(spline_factors,
         function(f) results(dds.jaw, name=paste0(f,".jaw",sp))[,2])
}
scjaw.md <- spline_coefs("md")
scjaw.mx <- spline_coefs("mx")
row.names(scjaw.md)=row.names(dds.jaw)
row.names(scjaw.mx)=row.names(dds.jaw)

scjaw.mus <- results(dds.jaw, name=paste0("species","mus"))[,2]
scjaw.ham <- results(dds.jaw, name=paste0("species","ham"))[,2]
names(scjaw.mus)=row.names(dds.jaw)
names(scjaw.ham)=row.names(dds.jaw)



## coeffs of the model with up/low in mouse, not in hamster. one intercept per species, and 3 knots per curve (upper mouse, lower mouse, hamster)

dds.3Mup <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethMusUp) 
dds.3Mup <- DESeq(dds.3Mup)

spline_coefs <- function(sp,type) {
  sapply(spline_factors,
         function(f) results(dds.3Mup, name=paste0(f,".paste.species..isupMus.",sp,".",type))[,2])
}
scjaw3Mup.musInt <- results(dds.3Mup, name=paste0("species","mus"))[,2]
scjaw3Mup.hamInt <- results(dds.3Mup, name=paste0("species","ham"))[,2]
names(scjaw3Mup.musInt)=row.names(dds.3Mup)
names(scjaw3Mup.hamInt)=row.names(dds.3Mup)

scjaw3Mup.mus.up <- spline_coefs("mus","yes")
scjaw3Mup.mus.lo <- spline_coefs("mus","no")
scjaw3Mup.ham <- spline_coefs("ham","no")
row.names(scjaw3Mup.mus.up)=row.names(dds.3Mup)
row.names(scjaw3Mup.mus.lo)=row.names(dds.3Mup)
row.names(scjaw3Mup.ham)=row.names(dds.3Mup)


## coeffs of the model with up/low in hamster, not in mouse
dds.3Hap <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethHamUp  ) 
dds.3Hap <- DESeq(dds.3Hap)

spline_coefs <- function(sp,type) {
  sapply(spline_factors,
         function(f) results(dds.3Hap, name=paste0(f,".paste.species..isupHam.",sp,".",type))[,2])
}
scjaw3Hap.musInt <- results(dds.3Hap, name=paste0("species","mus"))[,2]
scjaw3Hap.hamInt <- results(dds.3Hap, name=paste0("species","ham"))[,2]
names(scjaw3Hap.musInt)=row.names(dds.3Hap)
names(scjaw3Hap.hamInt)=row.names(dds.3Hap)

scjaw3Hap.ham.yes <- spline_coefs("ham","yes")
scjaw3Hap.ham.no <- spline_coefs("ham","no")
scjaw3Hap.mus <- spline_coefs("mus","no")
row.names(scjaw3Hap.ham.yes)=row.names(dds.3Hap)
row.names(scjaw3Hap.ham.no)=row.names(dds.3Hap)
row.names(scjaw3Hap.mus)=row.names(dds.3Hap)


## coeffs of the simple model

dds.simple <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.simple  ) 
dds.simple <- DESeq(dds.simple)

spline_coefs <- 
  sapply(spline_factors,
         function(f) results(dds.simple, name=f)[,2])

sc1teeth <- spline_coefs
row.names(sc1teeth)=row.names(dds.simple)

sc1teeth.mus <- results(dds.simple, name=paste0("species","mus"))[,2]
sc1teeth.ham <- results(dds.simple, name=paste0("species","ham"))[,2]
names(sc1teeth.mus)=row.names(dds.simple)
names(sc1teeth.ham)=row.names(dds.simple)

dds.diff.sptosimple <- DESeq(dds.diff.sp, test="LRT", reduced=design.simple)
res.diff.sptosimple <- results(dds.diff.sptosimple )

dds.jawtosimple <- DESeq(dds.jaw, test="LRT", reduced=design.simple)
res.jawtosimple <- results(dds.jawtosimple )

dds.sptosimple <- DESeq(dds.sp, test="LRT", reduced=design.simple)
res.sptosimple <- results(dds.sptosimple )



################################################################### plot full models in transcription factors

design.species.jaw.full <- as.formula(paste(" ~species*jaw  + (",
                                                paste(spline_factors,collapse=" + "),
                                                "):paste(species,jaw) + 0 "))
ddsfull <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.species.jaw.full ) 
ddsfull <- DESeq(ddsfull)

plotgenestage4cfull=function(id="Shh")
{
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  times = seq(-1,max(df1$timerel),length.out=100)
  k <- counts(ddsfull)[id,] / sizeFactors(ddsfull)
  dats=data.frame(times=df1$timerel,value=k,jaw=df1$jaw,species=df1$species)
  spline_coefs <- function(sp,mac) {
    sapply(spline_factors,
           function(f) results(ddsfull, name=paste0(f,".paste.species..jaw.",sp,".",mac))[,2])
  }
  sc.ham.md <- spline_coefs("ham","md")
  sc.ham.mx <- spline_coefs("ham","mx")
  sc.mus.md <- spline_coefs("mus","md")
  sc.mus.mx <- spline_coefs("mus","mx")
  row.names(sc.ham.md)=row.names(dds)
  row.names(sc.ham.mx)=row.names(dds)
  row.names(sc.mus.md)=row.names(dds)
  row.names(sc.mus.mx)=row.names(dds)
  
  
  sc.mus <- results(ddsfull, name=paste0("species","mus"))[,2]
  sc.ham <- results(ddsfull, name=paste0("species","ham"))[,2]
  names(sc.mus)=row.names(ddsfull)
  names(sc.ham)=row.names(ddsfull)
  sc.mx <- results(ddsfull, name=paste0("jaw","mx"))[,2]
  sc.mxmus <- results(ddsfull, name=paste0("speciesmus.jawmx"))[,2]
  names(sc.mx)=row.names(ddsfull)
  names(sc.mxmus)=row.names(ddsfull)
  
  
  
  # 4 different curves
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)
  dat$pred.ham.md <- as.numeric(2 ** (sc.ham[id]+sc.ham.md[id,] %*% SDM))
  dat$pred.ham.mx <- as.numeric(2 ** (sc.ham[id]+sc.mx[id]+sc.ham.mx[id,] %*% SDM))
  dat$pred.mus.md <- as.numeric(2 ** (sc.mus[id]+sc.mus.md[id,] %*% SDM))
  dat$pred.mus.mx <- as.numeric(2 ** (sc.mus[id]+sc.mx[id]+sc.mxmus[id]+sc.mus.mx[id,] %*% SDM))
  dat_melt=melt(dat,id.vars="times")
  dat_melt$jaw[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.mus.md"]="md"
  dat_melt$jaw[dat_melt$variable=="pred.ham.mx"|dat_melt$variable=="pred.mus.mx"]="mx"
  dat_melt$species[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.ham.mx"]="ham"
  dat_melt$species[dat_melt$variable=="pred.mus.md"|dat_melt$variable=="pred.mus.mx"]="mus"
  dat_melt$type="4 curves"
  
  
  models=data.frame(dat_melt)
  
  colo=c("ham md"="dark grey",
         "ham mx"="black",
         "mus md"="dark grey",
         "mus mx"="black")
  
  models=models[models$times>=0,]
  
  p2=ggplot(dats, aes(x=times, y=value, group= paste(species,jaw), col = paste(species,jaw))) + 
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
  #guides(color=FALSE)
  
  return(p2)
  
}



listOfGenes=read.csv("ListOfGenes.csv",h=F,sep=",")
system("mkdir ListOfGenes")

lapply(1:nrow(listOfGenes),function(r){
s=listOfGenes[r,]
id1=s[1,1]
if(id1%in%row.names(resnested))
{
a=plotgenestage4cfull(id=id1)+ theme(legend.position="none")
ggsave(paste0("ListOfGenes/",id1,"_",s[1,2],".pdf"),a,width = 8, height = 10, units = "cm") 
}
if(!id1%in%row.names(resnested))
{
  cat(paste0(id1,"\n"))
}
})

id1="Fgf10"
a=plotgenestage4cfull(id=id1)+ theme(legend.position="none")
ggsave(paste0("ListOfGenes/",id1,"_",s[1,2],".pdf"),a,width = 8, height = 10, units = "cm") 




listFT=c("Dlx1","Dlx5","Dlx6","Rarb","Barx1","Irx1","Lhx6","Lhx8",
"Mitf","Pou3f3","Tbx2","Lef1","Pitx1","Pax9","Six1", 
"Six2","Sox5","Sox8","Msx1","Osr2","Osr1","Foxd2",
"Hey2","Barx2","Alx1","Foxp2","Pbx1","Nkx3.1","Nkx2.3")

sapply(listFT,function(id1){
  if(id1%in%row.names(resnested))
  {
    a=plotgenestage4cfull(id=id1)+ theme(legend.position="none")
    ggsave(paste0("listFT/",id1,"_",s[1,2],".pdf"),a,width = 8, height = 10, units = "cm") 
  }
  if(!id1%in%row.names(resnested))
  {
    cat(paste0(id1,"\n"))
  }
})



####### to reproduce the supplementary figure with selected transcription factors

listFT=c("Nkx2.3","Pou3f3","Dlx1","Dlx5","Dlx6","Lhx6","Lhx8","Barx1","Pitx1",
         "Pax9","Msx1","Osr1","Osr2","Irx1","Six1","Hey2","Alx1","Rarb","Mitf","Lef1")
listFT=sort(listFT)         

l=lapply(sort(listFT),function(id1){
    a=plotgenestage4cfull(id=id1) + theme(legend.position="none")
    a
})
names(l)=listFT
p <- plot_grid(l, tags = listFT,margin = c(0.1, 0.1, 0.1, 0.1))
  ggsave("figexamples_forSupFig.png",p,width = 30, height = 30, units = "cm") 
  ggsave("figexamples_forSupFig.pdf",p,width = 40, height = 30, units = "cm") 
  



################################################################### study of models with 2 molars in one species

# compares full and up/low in mouse, not in hamster
  dds.4jaw3teethMusUp <- DESeq(dds.4teeth , test="LRT", reduced=design.3teethMusUp )
  res.diff.4jaw3teethMusUp <- results(dds.4jaw3teethMusUp)
  
# compares full and up/low in hamster, not in mouse
  dds.4jaw3teethHamUp <- DESeq(dds.4teeth , test="LRT", reduced=design.3teethHamUp )
  res.diff.4jaw3teethHamUp <- results(dds.4jaw3teethHamUp)
  
### summary of all results
  resnested=data.frame(c4vssp=res.diff.pat_species[,"padj"],
                       c4vsTooth=res.diff.pat_jaw[,"padj"], 
                       Spvs1=res.sptosimple[,"padj"],
                       Tovs1=res.jawtosimple[,"padj"],
                       Ham3t=res.diff.3teethHamUp[,"padj"],
                       Mus3t=res.diff.3teethMusUp[,"padj"],
                       Mus4t3t=res.diff.4jaw3teethMusUp[,"padj"],
                       Ham4t3t=res.diff.4jaw3teethHamUp[,"padj"]
  )
  row.names(resnested)=row.names(res.diff.pat_species)
  
  ## choice of the models: 
  resnested$type="pas"
  resnested$type[resnested$Spvs1<0.05&resnested$Tovs1>0.05]="Sp"
  resnested$type[resnested$Spvs1>0.05&resnested$Tovs1<0.05]="Tooth"
  resnested$type[resnested$Spvs1>0.05&resnested$Tovs1>0.05]="1 curve"
  resnested$type[resnested$type=="Sp"&resnested$c4vssp<0.05]="4 curves"
  resnested$type[resnested$type=="Tooth"&resnested$c4vsTooth<0.05]="4 curves"
  resnested$type[resnested$Spvs1<0.05&resnested$Tovs1<0.05&resnested$c4vsTooth<0.05&resnested$c4vssp<0.05]="4 curves"
  resnested$type[is.na(resnested$Spvs1)|is.na(resnested$Tovs1)|is.na(resnested$c4vssp)]="na"


###################################### classify genes that prefer 3 curves models over 2 species models
model3teeth=data.frame(cbind(Ham3teethpadj=res.diff.3teethHamUp$padj,Mus3teethpadj=res.diff.3teethMusUp$padj,Spvs1=res.sptosimple$padj,Tovs1=res.mactosimple$padj))
row.names(model3teeth)=row.names(res.diff.3teethHamUp)
model3teeth$model3teeth="1 curve"
model3teeth$model3teeth[res.diff.3teethHamUp$padj<0.05]="Ham3t"
model3teeth$model3teeth[res.diff.3teethMusUp$padj<0.05]="Mus3t"
model3teeth$model3teeth[res.diff.3teethMusUp$padj<0.05&res.diff.3teethHamUp$padj<0.05]="both3t"

###################################### rank genes in these models 
m=model3teeth[order(model3teeth$Mus3teethpadj),]
m$rank=1:nrow(model3teeth)
write.table(file="model3teeth.txt",model3teeth,sep="\t",row.names=T,col.names=T)

#+ Bucco/lingual differential expression for BMPER as cited in paper
#bl=read.table("../../DeconvolutionsBL/resultsBL28Juil21.txt",h=T)
#bl=bl[order(bl$padj),]
#bl$rank=1:nrow(bl)


###################################### PLOTS used in paper Figure 2

plotgenestage3tb=function(id="Bmper")
{
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  k <- counts(dds)[id,] / sizeFactors(dds)
  dats=data.frame(times=df1$timerel,value=k,jaw=df1$jaw,species=df1$species)
  SDM <- t(spline_design_matrix(times, intercept=F))
  
  # 3curves Ham
  dat3=data.frame(times=times)
  dat3$pred.ham.md <- as.numeric(2 ** (scjaw3Hap.hamInt[id]+scjaw3Hap.ham.no[id,] %*% SDM))
  dat3$pred.ham.mx <- as.numeric(2 ** (scjaw3Hap.hamInt[id]+scjaw3Hap.ham.yes[id,] %*% SDM))
  dat3$pred.mus.md <- as.numeric(2 ** (scjaw3Hap.musInt[id]+scjaw3Hap.mus[id,] %*% SDM))
  dat3$pred.mus.mx <- as.numeric(2 ** (scjaw3Hap.musInt[id]+scjaw3Hap.mus[id,] %*% SDM))
  dat3_melt=melt(dat3,id.vars="times")
  dat3_melt$jaw[dat3_melt$variable=="pred.ham.md"|dat3_melt$variable=="pred.mus.md"]="md"
  dat3_melt$jaw[dat3_melt$variable=="pred.ham.mx"|dat3_melt$variable=="pred.mus.mx"]="mx"
  dat3_melt$species[dat3_melt$variable=="pred.ham.md"|dat3_melt$variable=="pred.ham.mx"]="ham"
  dat3_melt$species[dat3_melt$variable=="pred.mus.md"|dat3_melt$variable=="pred.mus.mx"]="mus"
  dat3_melt$type="Sp + 2teeth ham"
  
  # 3curves Mus
  dat4=data.frame(times=times)
  dat4$pred.ham.md <- as.numeric(2 ** (scjaw3Mup.hamInt[id]+scjaw3Mup.ham[id,] %*% SDM))
  dat4$pred.ham.mx <- as.numeric(2 ** (scjaw3Mup.hamInt[id]+scjaw3Mup.ham[id,] %*% SDM))
  dat4$pred.mus.md <- as.numeric(2 ** (scjaw3Mup.musInt[id]+scjaw3Mup.mus.lo[id,] %*% SDM))
  dat4$pred.mus.mx <- as.numeric(2 ** (scjaw3Mup.musInt[id]+scjaw3Mup.mus.up[id,] %*% SDM))
  dat4_melt=melt(dat4,id.vars="times")
  dat4_melt$jaw[dat4_melt$variable=="pred.ham.md"|dat4_melt$variable=="pred.mus.md"]="md"
  dat4_melt$jaw[dat4_melt$variable=="pred.ham.mx"|dat4_melt$variable=="pred.mus.mx"]="mx"
  dat4_melt$species[dat4_melt$variable=="pred.ham.md"|dat4_melt$variable=="pred.ham.mx"]="ham"
  dat4_melt$species[dat4_melt$variable=="pred.mus.md"|dat4_melt$variable=="pred.mus.mx"]="mus"
  dat4_melt$type="Sp + 2teeth mus"

  # one curve per species
  dat2=data.frame(times=times)
  dat2$pred.ham.md <- as.numeric(2 ** (sc.sp.fham[id]+sc.sp.ham[id,] %*% SDM))
  dat2$pred.ham.mx <- as.numeric(2 ** (sc.sp.fham[id]+sc.sp.ham[id,] %*% SDM))
  dat2$pred.mus.md <- as.numeric(2 ** (sc.sp.fmus[id]+sc.sp.mus[id,] %*% SDM))
  dat2$pred.mus.mx <- as.numeric(2 ** (sc.sp.fmus[id]+sc.sp.mus[id,] %*% SDM))
  dat2_melt=melt(dat2,id.vars="times")
  dat2_melt$jaw[dat2_melt$variable=="pred.ham.md"|dat2_melt$variable=="pred.mus.md"]="md"
  dat2_melt$jaw[dat2_melt$variable=="pred.ham.mx"|dat2_melt$variable=="pred.mus.mx"]="mx"
  dat2_melt$species[dat2_melt$variable=="pred.ham.md"|dat2_melt$variable=="pred.ham.mx"]="ham"
  dat2_melt$species[dat2_melt$variable=="pred.mus.md"|dat2_melt$variable=="pred.mus.mx"]="mus"
  dat2_melt$type="species"
  
  models=data.frame(rbind(dat2_melt,dat3_melt,dat4_melt))
  
  
  p=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species ~ jaw, labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("species"="grey",
                                 "Sp + 2teeth mus"="purple",
                                 "Sp + 2teeth ham"="green"))+
    # guides(alpha=FALSE,color=FALSE)+
    ggtitle(paste0(id,"4vssp=",format(res.diff.pat_species[id,"padj"],digit=2),
                   " 4vsTooth=",format(res.diff.pat_jaw[id,"padj"],digit=2), 
                   " Spvs1=",format(res.sptosimple[id,"padj"],digit=2),
                   " 3tH=",format(res.diff.3teethHamUp[id,"padj"],digit=2),
                   " 3tM=",format(res.diff.3teethMusUp[id,"padj"],digit=2)))+
    theme(plot.title = element_text(size = 12, face = "bold"))
  cowplot::save_plot(paste0(id,"rel_3teeth.pdf"),p,
  ncol=1.5,nrow=2,
            limitsize = FALSE )     
  return(p)
}


a=plotgenestage3tb("Bmper")
b=plotgenestage3tb("Sfrp2")
c=plotgenestage3tb("Prickle1")
d=plotgenestage3tb("Cyp26b1")
e=plotgenestage3tb("Wnt11")


########################" sort genes into classes and plot for Figure 2C

getclassesmod3=function(list=biteit$V1,suf="bite-it",threshold=0.05)
{
  subset=model3teeth[row.names(model3teeth)%in%list,]
  out=c(table(subset$model),nrow(subset),suf)
  names(out)=c(names(table(subset$model)),"N","subset")
  out=out[c("subset","N","1 curve","both3t","Ham3t","Mus3t")]
  return(out)
}


n1=getclassesmod3(list=biteit$V1,suf="bite-it")
n2=getclassesmod3(list=dispensable$Gene..name,suf="Dispensable")
n3=getclassesmod3(list=keystone$Gene.name,suf="Keystone")
n6=getclassesmod3(list=row.names(model3teeth),suf="Total")
n7=getclassesmod3(list=pathwaylist,suf="Pathways")


datvar=data.frame(rbind(n1,n2,n3,n6,n7))
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

datvarmelt2$divham="no"
datvarmelt2$divham[datvarmelt2$variable=="Ham3t"|datvarmelt2$variable=="both3t"]="ham"

datvarmelt2$divmus="no"
datvarmelt2$divmus[datvarmelt2$variable=="Mus3t"|datvarmelt2$variable=="both3t"]="mus"



theme_barplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())



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


################## For FIgure 2C

propclass=ggplot(datvarmelt2species, aes(x = freq, y = subset3,fill = type)) +
  geom_bar(position=position_dodge(),stat="identity")+
  theme_barplot+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Ham 2curves","Mus 2curves"))

ggsave("SplineDivergentClassesBarplots.png",propclass,width = 12, height = 7, units = "cm") 
ggsave("SplineDivergentClassesBarplots.pdf",propclass,width = 12, height = 7, units = "cm")


### Proportion per signalling pathways

r=data.frame(lapply(names(pathwaylistspe),function(x){
  v=pathwaylistspe[[x]]
  n=getclassesmod3(list=v,suf=x,threshold=0.05)
  names(n)=c("subset","N","1 curve","both 3t","ham 2 teeth","mus 2 teeth")
  n[is.na(n)]=0
  return(n)
}))
names(r)=names(pathwaylistspe)
r
datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmeltp=melt(datvar,id=c("subset","N"))
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

## barplots per signalling pathway (not shown in paper)

datvarmelt2p$divham="no"
datvarmelt2p$divham[datvarmelt2p$variable=="Ham3t"|datvarmelt2p$variable=="both3t"]="ham"

datvarmelt2p$divmus="no"
datvarmelt2p$divmus[datvarmelt2p$variable=="Mus3t"|datvarmelt2p$variable=="both3t"]="mus"


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
  scale_fill_viridis(name="Model",discrete = TRUE, option = "C",labels = c("Ham 3curves","Mus 3curves"))

ggsave("SplineDivergent3ClassesPatBarplots.png",propclass,width = 15, height = 7, units = "cm") 
ggsave("SplineDivergent3ClassesPatBarplots.pdf",propclass,width = 15, height = 7, units = "cm")



################################################ GO analysis in supplementary figure 

## custom background = total 14532 gene list (with orthologs hamster/mouse)

lists=list(rownames(model3teeth)[model3teeth$model3teeth=="Mus3t"],rownames(model3teeth)[model3teeth$model3teeth=="Ham3t"])
names(lists)=c("mus","ham")

gostres <- gost(query = lists, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = rownames(model3teeth), 
                numeric_ns = "", sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"), as_short_link = FALSE)


#gostplot(gostres, capped = TRUE, interactive = TRUE)

gem <- gostres$result[,c("query","term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("model","GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem <- gem[,c("model","GO.ID", "Description", "FDR")]
write.table(gem, file = "gProfiler_gem_Mus3Ham3Models.txt", sep = "\t", quote = F, row.names = F)

set.seed(888)
go_id = gostres$result$term_id[gostres$result$query=="ham"]
mat = GO_similarity(go_id,ont="BP")
df = simplifyGO(mat)
pdf("simplifyGO_ham.pdf")
simplifyGO(mat)
dev.off()


go_id = gostres$result$term_id[gostres$result$query=="mus"]
mat = GO_similarity(go_id,ont="BP")
df = simplifyGO(mat)
pdf("simplifyGO_mus.pdf")
simplifyGO(mat)
dev.off()








