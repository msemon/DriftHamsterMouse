##### MODELS computed on the 4 teeth taken together. 

# load libraries
setwd("~/Documents/Projet_Drift/MappingGenome/Counts/Spline_Final/TimeGAM4Curves")
library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(readxl)
library(forcats)


## load whole tooth data
# metadata
metadataTot1=read.table(file="metadataTot.txt",sep="\t")
names(metadataTot1)=c("jaw","stage","species","rep","file")
# raw counts for 64 samples
CountTot=read.table(file="CountTot.txt",sep="\t")
CountTot=CountTot[,metadataTot1$jaw=="mx"|metadataTot1$jaw=="md"]
metadataTot1=metadataTot1[metadataTot1$jaw=="mx"|metadataTot1$jaw=="md",]

metadataTot1$stage[metadataTot1$stage==12.2]=12.25
metadataTot1$ech=paste(metadataTot1$species,metadataTot1$jaw,metadataTot1$stage,metadataTot1$rep,sep="_")
metadataTot1$id=row.names(metadataTot1)
metadataTot1=metadataTot1[,!names(metadataTot1)%in%c("stage")]


## time estimates from GAM model
mdhamtime=read.table("../predictions_rnaseq_cuspmdH_all.txt",h=T)
mxhamtime=read.table("../predictions_rnaseq_cuspmxH_all.txt",h=T)
mdmustime=read.table("../predictions_rnaseq_cuspmdM_all.txt",h=T)
mxmustime=read.table("../predictions_rnaseq_cuspmxM_all.txt",h=T)
times=data.frame(rbind(mdhamtime,mxhamtime,mdmustime,mxmustime))
names(times)=c("samples","stage","replicate","weight","fit_boxCox","lwr_boxCox","upr_boxCox","fit_log","lwr_log","upr_log","est_GAM","lower_GAM","upper_GAM")
times$sample[grep("mus",times$sample)]=paste0(times$sample[grep("mus",times$sample)],"W")

# add these estimated development times to colData and reorder samples in countdata
metadataTot2=merge(metadataTot1,times,by.x=0,by.y="sample",all.x=T,no.dups=TRUE)
CountTot2=CountTot[,metadataTot2$id]

ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot2,0),
                                        colData = metadataTot2,
                                        design = ~1)

colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="mus"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="ham"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="ham"])-12)/(14.5-12)
ddsToothWhole <- DESeq(ddsToothWhole)

### lists of genes of interest 
biteit=read.table("../liste_bite-it.csv",h=F)
keystone=read.csv("../keystone_genes.csv",h=T,sep="\t")
dispensable=read.csv("../dispensable_genes.csv",h=T,sep="\t")
pathway = read_excel("../pathways-Margaux-corMS.xlsx",sheet=1)
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
### make DESeq object
save(file="ddsToothWhole.Rdata",ddsToothWhole)


# simple plot gene by gene (no model fit)
plotgeneStage=function(id="Shh")
{
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  d <- plotCounts(ddsToothWhole, gene=id, intgroup=c("species","jaw","stage","timerel","rep"), 
                  returnData=TRUE)
   p=ggplot(d, aes(x=timerel, y=count, group=species:jaw,shape=rep)) + 
    facet_grid(species ~ jaw, labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2) + theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)") + guides(shape=FALSE)
  save_plot(paste0(id,"rel.png"), ncol=2,nrow=2,
            p ,limitsize = FALSE )    
}
plotgeneStage("Pou3f3")




###### Fit spline functions
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



## coeffs of the model with up/low in mouse, not in hamster

dds.3Mup <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethMusUp  ) 
dds.3Mup <- DESeq(dds.3Mup)

spline_coefs <- function(sp) {
  sapply(spline_factors,
         function(f) results(dds.3Mup, name=paste0(f,".paste.species..isupMus.mus.",sp))[,2])
}
scjaw3Mup.no <- spline_coefs("no")
scjaw3Mup.yes <- spline_coefs("yes")
row.names(scjaw3Mup.no)=row.names(dds.3Mup)
row.names(scjaw3Mup.yes)=row.names(dds.3Mup)

scjaw3Mup.mus <- results(dds.3Mup, name=paste0("species","mus"))[,2]
scjaw3Mup.ham <- results(dds.3Mup, name=paste0("species","ham"))[,2]
names(scjaw3Mup.mus)=row.names(dds.3Mup)
names(scjaw3Mup.ham)=row.names(dds.3Mup)


## coeffs of the model with up/low in hamster, not in mouse

dds.3Hap <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design = design.3teethHamUp  ) 
dds.3Hap <- DESeq(dds.3Hap)

spline_coefs <- function(sp) {
  sapply(spline_factors,
         function(f) results(dds.3Hap, name=paste0(f,".paste.species..isupHam.ham.",sp))[,2])
}
scjaw3Hap.no <- spline_coefs("no")
scjaw3Hap.yes <- spline_coefs("yes")
row.names(scjaw3Hap.no)=row.names(dds.3Hap)
row.names(scjaw3Hap.yes)=row.names(dds.3Hap)

scjaw3Hap.mus <- results(dds.3Hap, name=paste0("species","mus"))[,2]
scjaw3Hap.ham <- results(dds.3Hap, name=paste0("species","ham"))[,2]
names(scjaw3Hap.mus)=row.names(dds.3Hap)
names(scjaw3Hap.ham)=row.names(dds.3Hap)

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


## plot all models for a given gene

plotgenestage=function(id="Shh")
{
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  k <- counts(dds)[id,] / sizeFactors(dds)
  dats=data.frame(times=df1$timerel,value=k,jaw=df1$jaw,species=df1$species)
  
  
  # 4 different curves
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)
  dat$pred.ham.md <- as.numeric(2 ** (sc.ham[id]+sc.ham.md[id,] %*% SDM))
  dat$pred.ham.mx <- as.numeric(2 ** (sc.ham[id]+sc.ham.mx[id,] %*% SDM))
  dat$pred.mus.md <- as.numeric(2 ** (sc.mus[id]+sc.mus.md[id,] %*% SDM))
  dat$pred.mus.mx <- as.numeric(2 ** (sc.mus[id]+sc.mus.mx[id,] %*% SDM))
  dat_melt=melt(dat,id.vars="times")
  dat_melt$jaw[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.mus.md"]="md"
  dat_melt$jaw[dat_melt$variable=="pred.ham.mx"|dat_melt$variable=="pred.mus.mx"]="mx"
  dat_melt$species[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.ham.mx"]="ham"
  dat_melt$species[dat_melt$variable=="pred.mus.md"|dat_melt$variable=="pred.mus.mx"]="mus"
  dat_melt$type="4 teeth"
  
  
  # one curve per tooth type
  dat1=data.frame(times=times)
  dat1$pred.ham.md <- as.numeric(2 ** (scjaw.ham[id]+scjaw.md[id,] %*% SDM))
  dat1$pred.ham.mx <- as.numeric(2 ** (scjaw.ham[id]+scjaw.mx[id,] %*% SDM))
  dat1$pred.mus.md <- as.numeric(2 ** (scjaw.mus[id]+scjaw.md[id,] %*% SDM))
  dat1$pred.mus.mx <- as.numeric(2 ** (scjaw.mus[id]+scjaw.mx[id,] %*% SDM))
  dat1_melt=melt(dat1,id.vars="times")
  dat1_melt$jaw[dat1_melt$variable=="pred.ham.md"|dat1_melt$variable=="pred.mus.md"]="md"
  dat1_melt$jaw[dat1_melt$variable=="pred.ham.mx"|dat1_melt$variable=="pred.mus.mx"]="mx"
  dat1_melt$species[dat1_melt$variable=="pred.ham.md"|dat1_melt$variable=="pred.ham.mx"]="ham"
  dat1_melt$species[dat1_melt$variable=="pred.mus.md"|dat1_melt$variable=="pred.mus.mx"]="mus"
  dat1_melt$type="Upper-lower"
  
  
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
  dat2_melt$type="Hamster-mouse"
  
  
  # one curve 
  dat0=data.frame(times=times)
  dat0$pred.ham.md <- as.numeric(2 ** (sc1teeth.ham[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.ham.mx <- as.numeric(2 ** (sc1teeth.ham[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.mus.md <- as.numeric(2 ** (sc1teeth.mus[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.mus.mx <- as.numeric(2 ** (sc1teeth.mus[id]+sc1teeth[id,] %*% SDM))
  dat0_melt=melt(dat0,id.vars="times")
  dat0_melt$jaw[dat0_melt$variable=="pred.ham.md"|dat0_melt$variable=="pred.mus.md"]="md"
  dat0_melt$jaw[dat0_melt$variable=="pred.ham.mx"|dat0_melt$variable=="pred.mus.mx"]="mx"
  dat0_melt$species[dat0_melt$variable=="pred.ham.md"|dat0_melt$variable=="pred.ham.mx"]="ham"
  dat0_melt$species[dat0_melt$variable=="pred.mus.md"|dat0_melt$variable=="pred.mus.mx"]="mus"
  dat0_melt$type="Simple"
  
  models=data.frame(rbind(dat0_melt,dat1_melt,dat2_melt,dat_melt))
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  #title gives adjusted pvalues for 4 tests
  # 4vssp : full model versus species model
  # 4vsTooth : full model versus jaw model
  # Spvs1 : species model versus simple 1 curve model
  # Tovs1 : jaw model versus simple 1 curve model
  
  p=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species ~ jaw, labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("Simple"="black",
                                 "Hamster-mouse"="blue",
                                 "Upper-lower"="red",
                                 "4 teeth"="yellow"))+
    ggtitle(paste0(id," 4vssp=",format(res.diff.pat_species[id,"padj"],digit=2),
                   " 4vsTooth=",format(res.diff.pat_jaw[id,"padj"],digit=2), 
                   " Spvs1=",format(res.sptosimple[id,"padj"],digit=2),
                   " Tovs1=",format(res.jawtosimple[id,"padj"],digit=2)))+
    theme(plot.title = element_text(size = 12, face = "bold"))
  save_plot(paste0(id,"rel.pdf"), ncol=2,nrow=3,
            p ,limitsize = FALSE )     
  
  save_plot(paste0(id,"rel.png"), ncol=1.5,nrow=2,
            p ,limitsize = FALSE )     
  
  
  return(p)
  
}

#### For Figure 1D
plotgenestage("Fgf3")
plotgenestage("Bmp4")
plotgenestage("Fgf10")
plotgenestage("Fzd1")
plotgenestage("Wnt5a")
plotgenestage("Bmper")
plotgenestage("Shh")

plotgenestage("Sox3")
plotgenestage("Shh")

#### to plot nested models in different plots

plotgenestageNested=function(id="Shh")
{
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  times = seq(0,max(df1$timerel),length.out=100)
  
  k <- counts(dds)[id,] / sizeFactors(dds)
  dats=data.frame(times=df1$timerel,value=k,jaw=df1$jaw,species=df1$species)
  
  
  # 4 different curves
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)
  dat$pred.ham.md <- as.numeric(2 ** (sc.ham[id]+sc.ham.md[id,] %*% SDM))
  dat$pred.ham.mx <- as.numeric(2 ** (sc.ham[id]+sc.ham.mx[id,] %*% SDM))
  dat$pred.mus.md <- as.numeric(2 ** (sc.mus[id]+sc.mus.md[id,] %*% SDM))
  dat$pred.mus.mx <- as.numeric(2 ** (sc.mus[id]+sc.mus.mx[id,] %*% SDM))
  dat_melt=melt(dat,id.vars="times")
  dat_melt$jaw[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.mus.md"]="md"
  dat_melt$jaw[dat_melt$variable=="pred.ham.mx"|dat_melt$variable=="pred.mus.mx"]="mx"
  dat_melt$species[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.ham.mx"]="ham"
  dat_melt$species[dat_melt$variable=="pred.mus.md"|dat_melt$variable=="pred.mus.mx"]="mus"
  dat_melt$type="4 curves"
  
  
  # one curve per tooth type
  dat1=data.frame(times=times)
  dat1$pred.ham.md <- as.numeric(2 ** (scjaw.ham[id]+scjaw.md[id,] %*% SDM))
  dat1$pred.ham.mx <- as.numeric(2 ** (scjaw.ham[id]+scjaw.mx[id,] %*% SDM))
  dat1$pred.mus.md <- as.numeric(2 ** (scjaw.mus[id]+scjaw.md[id,] %*% SDM))
  dat1$pred.mus.mx <- as.numeric(2 ** (scjaw.mus[id]+scjaw.mx[id,] %*% SDM))
  dat1_melt=melt(dat1,id.vars="times")
  dat1_melt$jaw[dat1_melt$variable=="pred.ham.md"|dat1_melt$variable=="pred.mus.md"]="md"
  dat1_melt$jaw[dat1_melt$variable=="pred.ham.mx"|dat1_melt$variable=="pred.mus.mx"]="mx"
  dat1_melt$species[dat1_melt$variable=="pred.ham.md"|dat1_melt$variable=="pred.ham.mx"]="ham"
  dat1_melt$species[dat1_melt$variable=="pred.mus.md"|dat1_melt$variable=="pred.mus.mx"]="mus"
  dat1_melt$type="Tooth"
  
  
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
  dat2_melt$type="Species"
  
  
  # one curve 
  dat0=data.frame(times=times)
  dat0$pred.ham.md <- as.numeric(2 ** (sc1teeth.ham[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.ham.mx <- as.numeric(2 ** (sc1teeth.ham[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.mus.md <- as.numeric(2 ** (sc1teeth.mus[id]+sc1teeth[id,] %*% SDM))
  dat0$pred.mus.mx <- as.numeric(2 ** (sc1teeth.mus[id]+sc1teeth[id,] %*% SDM))
  dat0_melt=melt(dat0,id.vars="times")
  dat0_melt$jaw[dat0_melt$variable=="pred.ham.md"|dat0_melt$variable=="pred.mus.md"]="md"
  dat0_melt$jaw[dat0_melt$variable=="pred.ham.mx"|dat0_melt$variable=="pred.mus.mx"]="mx"
  dat0_melt$species[dat0_melt$variable=="pred.ham.md"|dat0_melt$variable=="pred.ham.mx"]="ham"
  dat0_melt$species[dat0_melt$variable=="pred.mus.md"|dat0_melt$variable=="pred.mus.mx"]="mus"
  dat0_melt$type="1 curve"
  
  models=data.frame(rbind(dat0_melt,dat1_melt,dat2_melt,dat_melt))
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  
  p0=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species~ jaw , labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")
  
  p1=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species~ jaw , labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models[models$type=="1 curve",],aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "Species"="blue",
                                 "Tooth"="red",
                                 "4 curves"="yellow"),
                        labels=c("1 curve"="1 curve",
                                 "Species"="2 species",
                                 "Tooth"="2 teeth",
                                 "4 curves"="4 curves"))+guides(color=FALSE)
  
  
  p2=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species~ jaw , labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models[models$type=="1 curve"|models$type=="Tooth",],aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "Species"="blue",
                                 "Tooth"="red",
                                 "4 curves"="yellow"),
                        labels=c("1 curve"="1 curve",
                                 "Species"="2 species",
                                 "Tooth"="2 teeth",
                                 "4 curves"="4 curves"))+guides(color=FALSE)
  
  p3=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species~ jaw , labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models[models$type=="1 curve"|models$type=="Species",],aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "Species"="blue",
                                 "Tooth"="red",
                                 "4 curves"="yellow"),
                        labels=c("1 curve"="1 curve",
                                 "Species"="2 species",
                                 "Tooth"="2 teeth",
                                 "4 curves"="4 curves"))+guides(color=FALSE)
  
  p4=ggplot(dats, aes(x=times, y=value, group=species)) + 
    facet_grid(species~ jaw , labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models[models$type=="4 curves"|models$type=="Species",],aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "Species"="blue",
                                 "Tooth"="red",
                                 "4 curves"="green"),
                        labels=c("1 curve"="1 curve",
                                 "Species"="2 species",
                                 "Tooth"="2 teeth",
                                 "4 curves"="4 curves"))+guides(color=FALSE)
  save_plot(paste0(id,"_model_without_curve.png"), ncol=0.8,nrow=1,
            p0 )     
  save_plot(paste0(id,"_model_simple_1curve.pdf"), ncol=0.8,nrow=1,
            p1 )     
  save_plot(paste0(id,"_model_jaw_2teeth.png"), ncol=0.8,nrow=1,
            p2 )     
  save_plot(paste0(id,"_model_species_2species.png"), ncol=0.8,nrow=1,
            p3 )     
  save_plot(paste0(id,"_model_complex_4teeth.pdf"), ncol=0.8,nrow=1,
            p4 )   
}

plotgenestageNested("Bmp4")


#### To plot examples, superimposed

plotgenestage4c=function(id="Shh")
{
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  times = seq(0,max(df1$timerel),length.out=100)
  k <- counts(dds)[id,] / sizeFactors(dds)
  dats=data.frame(times=df1$timerel,value=k,jaw=df1$jaw,species=df1$species)
  
  # 4 different curves
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)
  dat$pred.ham.md <- as.numeric(2 ** (sc.ham[id]+sc.ham.md[id,] %*% SDM))
  dat$pred.ham.mx <- as.numeric(2 ** (sc.ham[id]+sc.ham.mx[id,] %*% SDM))
  dat$pred.mus.md <- as.numeric(2 ** (sc.mus[id]+sc.mus.md[id,] %*% SDM))
  dat$pred.mus.mx <- as.numeric(2 ** (sc.mus[id]+sc.mus.mx[id,] %*% SDM))
  dat_melt=melt(dat,id.vars="times")
  dat_melt$jaw[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.mus.md"]="md"
  dat_melt$jaw[dat_melt$variable=="pred.ham.mx"|dat_melt$variable=="pred.mus.mx"]="mx"
  dat_melt$species[dat_melt$variable=="pred.ham.md"|dat_melt$variable=="pred.ham.mx"]="ham"
  dat_melt$species[dat_melt$variable=="pred.mus.md"|dat_melt$variable=="pred.mus.mx"]="mus"
  dat_melt$type="4 curves"
  
  
  models=data.frame(dat_melt)
  
  colo=c("ham md"="light green",
         "ham mx"="dark green",
         "mus md"="light blue",
         "mus mx"="dark blue")
  colo=c("ham md"="green",
         "ham mx"="yellow",
         "mus md"="blue",
         "mus mx"="purple")
  
  p2=ggplot(dats, aes(x=times*10, y=value, group= paste(species,jaw), col = paste(species,jaw))) + 
    facet_grid(species ~., labeller = labeller(species=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=2,alpha = 0.2) + theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models,aes(x=times*10, y=value,group= paste(species,jaw),col = paste(species,jaw)))+
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

a=plotgenestage4c("Bmper")
b=plotgenestage4c("Fgf3")
c=plotgenestage4c("Fgf10")
d=plotgenestage4c("Bmp4")
e=plotgenestage4c("Wif1")
f=plotgenestage4c("Sox2")

ggsave("FGF10.pdf",c,width = 12, height = 8, units = "cm") 
ggsave("Bmper.pdf",a,width = 12, height = 8, units = "cm") 
ggsave("Fgf3.pdf",b,width = 12, height = 8, units = "cm") 
ggsave("Bmp4_a.pdf",d,width = 12, height = 8, units = "cm") 



prow2 <- plot_grid(
  b+ theme(legend.position="none"), 
  d+ theme(legend.position="none"),
  c+ theme(legend.position="none"),
  a+ theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  ncol = 1
)
legend2 <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 4))
)
figexamples <- plot_grid(prow2, legend2, rel_widths = c(3, 1))
ggsave("figexamples.png",figexamples,width = 8, height = 29, units = "cm") 






prow2 <- plot_grid(
  a+ theme(legend.position="none"), 
  e+ theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  ncol = 1
)
legend2 <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 4))
)
figexamples <- plot_grid(prow2, legend2, rel_widths = c(3, .7))
ggsave("figexamples2.png",figexamples,width = 8, height = 29, units = "cm") 





prow2 <- plot_grid(
  a+ theme(legend.position="none"), 
  d+ theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  ncol = 1
)
legend2 <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 4))
)
figexamples <- plot_grid(prow2, legend2, rel_widths = c(3, .7))
ggsave("figexamples3.png",figexamples,width = 8, height = 29, units = "cm") 



## get 3 teeth models by comparison to full model
# full model needs to be specified with factors in simpler model
design.4teethMusHam <- as.formula(paste("~ species + (",
                                        paste(spline_factors,collapse=" + "),
                                        "):paste(species,isupMus,isupHam) + 0 "))
dds.4teeth <- DESeqDataSetFromMatrix(countData =  counts(ddsToothWhole)[,colData(ddsToothWhole)$stage>=12], colData =df1 , design =  design.4teethMusHam) 

# compares full and up/low in mouse, not in hamster
dds.4teeth3teethMusUp <- DESeq(dds.4teeth , test="LRT", reduced=design.3teethMusUp )
res.diff.4teeth3teethMusUp <- results(dds.4teeth3teethMusUp)

# compares full and up/low in hamster, not in mouse
dds.4teeth3teethHamUp <- DESeq(dds.4teeth , test="LRT", reduced=design.3teethHamUp )
res.diff.4teeth3teethHamUp <- results(dds.4teeth3teethHamUp)



### summary of all results

resnested=data.frame(c4vssp=res.diff.pat_species[,"padj"],
                     c4vsTooth=res.diff.pat_jaw[,"padj"], 
                     Spvs1=res.sptosimple[,"padj"],
                     Tovs1=res.jawtosimple[,"padj"],
                     Ham3t=res.diff.3teethHamUp[,"padj"],
                     Mus3t=res.diff.3teethMusUp[,"padj"],
                     Mus4t3t=res.diff.4teeth3teethMusUp[,"padj"],
                     Ham4t3t=res.diff.4teeth3teethHamUp[,"padj"]
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


# compute the distance between teeth, using complete (4 curves) model (not used in paper but corroborates the BCA)
times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
SDM <- t(spline_design_matrix(times, intercept=F))

dat=data.frame(times=times)
dat0=data.frame(times=times)

ldc=sapply(row.names(resnested),function(gene){
  
  dat=data.frame(times=times)
  dat$pred.ham.md <- as.numeric(2 ** (sc.ham[gene]+sc.ham.md[gene,] %*% SDM))
  dat$pred.ham.mx <- as.numeric(2 ** (sc.ham[gene]+sc.ham.mx[gene,] %*% SDM))
  dat$pred.mus.md <- as.numeric(2 ** (sc.mus[gene]+sc.mus.md[gene,] %*% SDM))
  dat$pred.mus.mx <- as.numeric(2 ** (sc.mus[gene]+sc.mus.mx[gene,] %*% SDM))
  
  # distance 
  dc=c(sqrt(sum((dat[,"pred.ham.md"]-dat[,"pred.ham.mx"])^2)),
       sqrt(sum((dat[,"pred.mus.md"]-dat[,"pred.mus.mx"])^2)),
       sqrt(sum((dat[,"pred.mus.md"]-dat[,"pred.ham.md"])^2)),
       sqrt(sum((dat[,"pred.mus.mx"]-dat[,"pred.ham.mx"])^2)))
  names(dc)=c("ham","mus","md","mx")
  
  # signe 
  dcS=c(sum(dat[,"pred.ham.md"]-dat[,"pred.ham.mx"])>0,
        sum(dat[,"pred.mus.md"]-dat[,"pred.mus.mx"])>0,
        sum(dat[,"pred.mus.md"]-dat[,"pred.ham.md"])>0,
        sum(dat[,"pred.mus.mx"]-dat[,"pred.ham.mx"])>0)
  names(dcS)=c("hamS","musS","mdS","mxS")
  
  
  # distance norm
  dc2=c(sqrt(sum((dat[,"pred.ham.md"]-dat[,"pred.ham.mx"])^2))/sqrt(sum((dat[,"pred.ham.md"]+dat[,"pred.ham.mx"])^2)),
        sqrt(sum((dat[,"pred.mus.md"]-dat[,"pred.mus.mx"])^2))/sqrt(sum((dat[,"pred.mus.md"]+dat[,"pred.mus.mx"])^2)),
        sqrt(sum((dat[,"pred.mus.md"]-dat[,"pred.ham.md"])^2))/sqrt(sum((dat[,"pred.mus.md"]+dat[,"pred.ham.md"])^2)),
        sqrt(sum((dat[,"pred.mus.mx"]-dat[,"pred.ham.mx"])^2))/sqrt(sum((dat[,"pred.mus.mx"]+dat[,"pred.ham.mx"])^2)))
  names(dc2)=c("hamN","musN","mdN","mxN")
  
  return(c(dc,dc2,dcS))
})

distances=data.frame(t(data.frame(ldc)))

par(mfrow=c(2,1))
pdf("distancesBmper_norm.pdf")
plot(distances$musN,distances$hamN)
points(distances["Bmper","musN"],distances["Bmper","hamN"],col="red",pch=20)
plot(distances$mdN,distances$mxN)
points(distances["Bmper","mdN"],distances["Bmper","mxN"],col="red",pch=20)
dev.off()

par(mfrow=c(2,1))
pdf("distancesBmper.pdf")
plot(distances$mus,distances$ham)
points(distances["Bmper","mus"],distances["Bmper","ham"],col="red",pch=20)
plot(distances$mdN,distances$mxN)
points(distances["Bmper","md"],distances["Bmper","mx"],col="red",pch=20)
dev.off()

resnested=cbind(resnested,distances)
write.table(resnested,file="resnested_TimeGAM.txt",row.names=T,sep="\t",col.names=T)


### Proportion of the 4 models per class. not used in paper
theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())

getclasses=function(list=biteit$V1,suf="bite-it",threshold=0.05)
{
  resnested$type="pas"
  resnested$type[resnested$Spvs1<threshold&resnested$Tovs1>threshold]="Sp"
  resnested$type[resnested$Spvs1>threshold&resnested$Tovs1<threshold]="Tooth"
  resnested$type[resnested$Spvs1>threshold&resnested$Tovs1>threshold]="1 curve"
  resnested$type[resnested$type=="Sp"&resnested$c4vssp<threshold]="4 curves"
  resnested$type[resnested$type=="Tooth"&resnested$c4vsTooth<threshold]="4 curves"
  resnested$type[resnested$Spvs1<threshold&resnested$Tovs1<threshold&resnested$c4vsTooth<threshold&resnested$c4vssp<threshold]="4 curves"
  resnested$type[is.na(resnested$Spvs1)|is.na(resnested$Tovs1)|is.na(resnested$c4vssp)]="na"
 
  subset=resnested[row.names(resnested)%in%list,]
  out=c(table(subset$type),nrow(subset),suf)
  names(out)=c(names(table(subset$type)),"N","subset")
  out=out[c("subset","N","1 curve","Sp","Tooth","4 curves")]
  return(out)
}

n1=getclasses(list=biteit$V1,suf="bite-it")
n2=getclasses(list=dispensable$Gene..name,suf="Dispensable")
n3=getclasses(list=keystone$Gene.name,suf="Keystone")
n6=getclasses(list=row.names(resnested),suf="Total")
n7=getclasses(list=pathwaylist,suf="Pathways")
datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
datvarmelt$N=as.numeric(as.character(datvarmelt$N))

datvarmelt$subset=factor(datvarmelt$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))

datvarmelt2 <- datvarmelt %>%
  group_by(subset) %>%
  mutate(percent = value / sum(value)) %>%
  ungroup()
datvarmelt2$subset=factor(datvarmelt2$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))


### Proportion per signalling pathways

r=data.frame(lapply(names(pathwaylistspe),function(x){
  v=pathwaylistspe[[x]]
  n=getclasses(list=v,suf=x,threshold=0.05)
  names(n)=c("subset","N","1 curve","Sp","Tooth","4 curves")
  n[is.na(n)]=0
  return(n)
}))
names(r)=names(pathwaylistspe)

datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
row.names(datvar)=datvar$subset
nvar=datvar$N
names(nvar)=datvar$subset

datvarp=data.frame(t(r))
datvarp=datvarp[order(as.numeric(datvarp$N),decreasing=TRUE),]
pator=datvarp$N
names(pator)=datvarp$subset


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
  scale_colour_manual(name="Model",
                      values=c("X1.curve"="grey",
                               "Sp"="blue",
                               "Tooth"="red",
                               "X4.curves"="purple"),
                      labels = c("1 curve", "Species", "Tooth","4 curves"))

propclass2 = propclass2 +guides(size=FALSE)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 13.2, ymin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 11.2, ymin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2$percent)+0.1,  ymax = 9.2, ymin = 0.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", xmin = 0, xmax = max(datvarmelt2p$percent)+0.1,  ymax = 9.2, ymin = 8.8, fill = "lightsteelblue", alpha = 0.2)
ggsave("SplineClassesPathways.png",propclass2,width = 15, height = 10, units = "cm") 

ggsave("SplineClassesPathways.pdf",propclass2,width = 15, height = 10, units = "cm") 


# with barplots

propclass2=ggplot(datvarmelt2p, aes(x = subset3, y=percent,  fill=variable)) +
  geom_bar(stat="identity", width=0.6,position=position_dodge())+
  theme_dotplot +
  scale_y_continuous(labels = scales::percent)+theme(axis.text.x = element_text(angle = 45))+coord_flip()+
  ylab("Gene subsets")+
  scale_fill_manual(name="Model",
                      values=c("X1.curve"="grey",
                               "Sp"="blue",
                               "Tooth"="red",
                               "X4.curves"="purple"),
                      labels = c("1 curve", "Species", "Tooth","4 curves"))

propclass2 = propclass2 + theme(legend.position = "none")+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 13.2, xmin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 11.2, xmin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 9.2, xmin = 0.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 9.2, xmin = 8.8, fill = "lightsteelblue", alpha = 0.2)

ggsave("SplineClassesPathways_bars.png",propclass2,width = 15, height = 10, units = "cm") 

ggsave("SplineClassesPathways_bars.pdf",propclass2,width = 15, height = 10, units = "cm") 





## Proportion coevolution measured as genes with at least 2 curves (even if they also have 2 teeth)
## In the paper corresponds to panel Figure 3C
## coevolution : genes better fit by species model than simple model.

resnested$coev="no"
resnested$coev[resnested$Spvs1<0.05]="Coev"

getclasses2=function(list=biteit$V1,suf="bite-it",threshold=0.05)
{
  resnested$type="no"
  resnested$type[resnested$Spvs1<threshold&resnested$Tovs1>threshold]="Sp"
  resnested$type[resnested$Spvs1>threshold&resnested$Tovs1<threshold]="Tooth"
  resnested$type[resnested$Spvs1>threshold&resnested$Tovs1>threshold]="1 curve"
  resnested$type[resnested$type=="Sp"&resnested$c4vssp<threshold]="4 curves"
  resnested$type[resnested$type=="Tooth"&resnested$c4vsTooth<threshold]="4 curves"
  resnested$type[resnested$Spvs1<threshold&resnested$Tovs1<threshold&resnested$c4vsTooth<threshold&resnested$c4vssp<threshold]="4 curves"
  resnested$type[is.na(resnested$Spvs1)|is.na(resnested$Tovs1)|is.na(resnested$c4vssp)]="na"
  
  subset=resnested[row.names(resnested)%in%list,]
  subset=subset[subset$type!="1 curve"&subset$type!="na",]
  out=c(table(subset$coev),nrow(subset),suf)
  names(out)=c(names(table(subset$coev)),"N","subset")
  out=out[c("subset","N","Coev","no")]
  return(out)
}

n1=getclasses2(list=biteit$V1,suf="bite-it")
n2=getclasses2(list=dispensable$Gene..name,suf="Dispensable")
n3=getclasses2(list=keystone$Gene.name,suf="Keystone")
n6=getclasses2(list=row.names(resnested),suf="Total")
n7=getclasses2(list=pathwaylist,suf="Pathways")
datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
datvarmelt$N=as.numeric(as.character(datvarmelt$N))
datvarmelt$subset=factor(datvarmelt$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))

datvarmelt2 <- datvarmelt %>%
  group_by(subset) %>%
  mutate(percent = value / sum(value))
datvarmelt2$subset=factor(datvarmelt2$subset,levels=c("Pathways","Keystone","Dispensable","bite-it","Total"))



r=data.frame(lapply(names(pathwaylistspe),function(x){
  v=pathwaylistspe[[x]]
  n=getclasses2(list=v,suf=x,threshold=0.05)
  names(n)=c("subset","N","Coev","no")
  n[is.na(n)]=0
  return(n)
}))

names(r)=names(pathwaylistspe)
datvar=data.frame(rbind(n1,n2,n3,n6,n7))
datvarmelt=melt(datvar,id=c("subset","N"))
datvarmelt$value=as.numeric(as.character(datvarmelt$value))
row.names(datvar)=datvar$subset
datvarp=data.frame(t(r))
datvarp=rbind(datvar[datvar$subset%in%c("bite-it","Pathways","Dispensable","Keystone","Total"),],datvarp)
datvarmeltp=melt(datvarp,id=c("subset","N"))
datvarmeltp$value=as.numeric(as.character(datvarmeltp$value))

datvarmelt2p <- datvarmeltp %>%
  group_by(subset) %>%
  mutate(percent = value / sum(value)) %>%
  ungroup()%>%
  filter(variable=="Coev")%>%
  mutate(subset2 = paste0(subset," (",N,")"))%>%
  mutate(subset2 = fct_reorder(subset2, as.numeric(as.character(N))))


propclass2=ggplot(datvarmelt2p, aes(x = subset2, y=percent)) +
  geom_bar(stat="identity", width=0.5)+
  theme_dotplot +
  scale_y_continuous(labels = scales::percent)+theme(axis.text.x = element_text(angle = 45))+coord_flip()

propclass2 = propclass2 + theme(legend.position = "none")+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 13.2, xmin = 12.8, fill = "grey", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 11.2, xmin = 9.8, fill = "palegreen", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 9.2, xmin = 0.8, fill = "lightsteelblue", alpha = 0.2)+
  annotate("rect", ymin = 0, ymax = max(datvarmelt2p$percent)+0.1,  xmax = 9.2, xmin = 8.8, fill = "lightsteelblue", alpha = 0.2)
ggsave("SplineCoevoPathways2.png",propclass2,width = 15, height = 10, units = "cm") 
ggsave("SplineCoevoPathways2.pdf",propclass2,width = 15, height = 10, units = "cm") 


