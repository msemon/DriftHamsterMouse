##### MODELS computed on pairs of teeth 

################################################################### load libraries

setwd("~/Documents/Projet_Drift/MappingGenome/Counts/Spline_Final/TimeGAM")
library(DESeq2)
library(splines)
library(tximport)
library(ggplot2)
library(ade4)
library(cowplot)
library(reshape2)
library(plyr)
library(readxl)


################################################################### load data
## load whole tooth data
# metadata
metadataTot1=read.table(file="../../metadataTot.txt",sep="\t")
names(metadataTot1)=c("jaw","stage","species","rep","file")
# raw counts for 64 samples
CountTot=read.table(file="../../CountTot.txt",sep="\t")
CountTot=CountTot[,metadataTot1$jaw=="mx"|metadataTot1$jaw=="md"]
metadataTot1=metadataTot1[metadataTot1$jaw=="mx"|metadataTot1$jaw=="md",]

metadataTot1$stage[metadataTot1$stage==12.2]=12.25
metadataTot1$ech=paste(metadataTot1$species,metadataTot1$jaw,metadataTot1$stage,metadataTot1$rep,sep="_")
metadataTot1$id=row.names(metadataTot1)
metadataTot1=metadataTot1[,!names(metadataTot1)%in%c("stage")]


## time estimates
mdhamtime=read.table("../predictions_rnaseq_cuspmdH_all.txt",h=T)
mxhamtime=read.table("../predictions_rnaseq_cuspmxH_all.txt",h=T)
mdmustime=read.table("../predictions_rnaseq_cuspmdM_all.txt",h=T)
mxmustime=read.table("../predictions_rnaseq_cuspmxM_all.txt",h=T)
times=data.frame(rbind(mdhamtime,mxhamtime,mdmustime,mxmustime))
names(times)=c("samples","stage","replicate","weight","fit_boxCox","lwr_boxCox","upr_boxCox","fit_log","lwr_log","upr_log","est_GAM","lower_GAM","upper_GAM")
times$sample[grep("mus",times$sample)]=paste0(times$sample[grep("mus",times$sample)],"W")

metadataTot2=merge(metadataTot1,times,by.x=0,by.y="sample",all.x=T,no.dups=TRUE)
CountTot2=CountTot[,metadataTot2$id]



### make DESeq object; we remove 12 young hamster samples
ddsToothWhole <- DESeqDataSetFromMatrix(countData = round(CountTot2,0)[,metadataTot2$stage>=12],
                              colData = metadataTot2[metadataTot2$stage>=12,],
                              design = ~1)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="mus"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="mus"])-14.5)/(18-14.5)
colData(ddsToothWhole)$timerel[colData(ddsToothWhole)$species=="ham"]=10*(as.numeric(colData(ddsToothWhole)$stage[colData(ddsToothWhole)$species=="ham"])-12)/(14.5-12)
ddsToothWhole <- DESeq(ddsToothWhole)





### lists of genes
biteit=read.table("../liste_bite-it.csv",h=F)
keystone=read.csv("../keystone_genes.csv",h=T,sep="\t")
dispensable=read.csv("../dispensable_genes.csv",h=T,sep="\t")
pathway = read_excel("../pathways-Margaux-corMS.xlsx",sheet=1)
pathway=pathway[,1:5]
pathwaylist=unique(pathway$Symbol)
pathwaylistspe=lapply(split(pathway,pathway$Pathway),function(x){unique(x$Symbol)})

##################################" compute the splines on the 2 selected timeseries
# if facselec="jaw" then computes species effect on one tooth (upper or lower)
# if facselec="species" then computes tooth effect on one species (mus or ham)

compute_timeeffect=function(facselec="jaw",selec="md",fac="species")
{
  dds=NULL
   res.diff=NULL
   sc.sp1=NULL
   sc.sp2=NULL
   sc.fac1=NULL
   sc.fac2=NULL
   sc.red.sp1=NULL
   sc.red.sp2=NULL
   sc.red.t=NULL
     
     
  dds <- DESeqDataSetFromMatrix(countData = counts(ddsToothWhole)[,colData(ddsToothWhole)[,facselec]==selec],
                                colData = colData(ddsToothWhole)[colData(ddsToothWhole)[,facselec]==selec,],
                                design = ~1)
  #### HERE we decide to take time estimate by GAM model
  colData(dds)$timerel[colData(dds)$species=="mus"]=10*(as.numeric(colData(dds)$est_GAM[colData(dds)$species=="mus"])-14.5)/(18-14.5)
  colData(dds)$timerel[colData(dds)$species=="ham"]=10*(as.numeric(colData(dds)$est_GAM[colData(dds)$species=="ham"])-12.3)/(14.5-12.3)
  dds <- DESeq(dds)
  
  spline_design_matrix <- function(time, intercept = F) {
    tmin <- min(time)
    tmax <- max(time)
    X <- bs(time,
            Boundary.knots = c(tmin,tmax),
            intercept = intercept,
            degree = 3)
    colnames(X) <- paste("bs_", colnames(X), sep = "")
    X
  }
  
  df1=colData(dds)
  SDM <- spline_design_matrix(df1$timerel, intercept = F)
  spline_factors <- colnames(SDM)
  
  if(facselec=="species")
  {
    design <- as.formula(paste(" ~ jaw +"," (",
                               paste(spline_factors,collapse=" + "),
                               ") : ",fac," + 0 "))
    
    red.design <- as.formula(paste(" ~","+"," (",
                                   paste( spline_factors,collapse=" + "),
                                   ") "))
   }
  if(facselec=="jaw")
  {  

  design <- as.formula(paste(" ~ species + (",
                             paste(spline_factors,collapse=" + "),
                             ") : ",fac," + 0 "))
  red.design <- as.formula(paste(" ~ species  + (",
                                 paste(spline_factors,collapse=" + "),
                                 ")  + 0 "))
  }
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts(dds),
                                colData = cbind(df1[,c("species","jaw","stage","rep","timerel")],SDM),
                                design = design)
  
  dds <- DESeq(dds)
  
  if(facselec=="species")
  {
  spline_coefs <- function(sp) {
    sapply(spline_factors,
           function(f) results(dds, name=paste0(fac,sp,".",f))[,2])
  }
  
  levelsfac=levels(as.factor(df1[,fac]))
  sc.fac1 <- spline_coefs(levelsfac[1])
  sc.fac2 <- spline_coefs(levelsfac[2])
  row.names(sc.fac1)=row.names(dds)
  row.names(sc.fac2)=row.names(dds)
  
  sc.sp1 <- results(dds, name=paste0("jaw","md"))[,2]
  sc.sp2 <- results(dds, name=paste0("jaw","mx"))[,2]
  names(sc.sp1)=row.names(dds)
  names(sc.sp2)=row.names(dds)
  }
  
  if(facselec=="jaw")
  {
    spline_coefs <- function(sp) {
      sapply(spline_factors,
             function(f) results(dds, name=paste0(fac,sp,".",f))[,2])
    }
    
    levelsfac=levels(as.factor(df1[,fac]))
    sc.fac1 <- spline_coefs(levelsfac[1])
    sc.fac2 <- spline_coefs(levelsfac[2])
    row.names(sc.fac1)=row.names(dds)
    row.names(sc.fac2)=row.names(dds)  
    
    
  sc.sp1 <- results(dds, name=paste0("species","ham"))[,2]
  sc.sp2 <- results(dds, name=paste0("species","mus"))[,2]
  names(sc.sp1)=row.names(dds)
  names(sc.sp2)=row.names(dds)
  }
  
  red.dds <- DESeqDataSetFromMatrix(countData=counts(dds),
                                    colData=cbind(df1[,c("jaw","stage","species","jaw","timerel")],SDM),
                                    design=red.design)
  
  red.dds <- DESeq(red.dds, betaPrior=F)
  sc.red.t <- sapply(spline_factors,function(f) results(red.dds, name=f)[,2])
  row.names(sc.red.t)=row.names(dds)
  
  if(facselec=="jaw")
  {
  sc.red.sp1 <- results(red.dds, name=paste0("species","ham"))[,2]
  sc.red.sp2 <- results(red.dds, name=paste0("species","mus"))[,2]
  names(sc.red.sp1)=row.names(red.dds)
  names(sc.red.sp2)=row.names(red.dds)
  }
  
  
   if(facselec=="species")
   {
     sc.red.sp1 <- results(red.dds, name="Intercept")[,2]
     sc.red.sp2 <- results(red.dds, name="Intercept")[,2]
     names(sc.red.sp1)=row.names(red.dds)
     names(sc.red.sp2)=row.names(red.dds)
   }
  
  dds.diff <- DESeq(dds, test="LRT", reduced=red.design)
  res.diff <- results(dds.diff)
  gene_ranking <- (1:nrow(res.diff))[order(res.diff$pvalue)]
  ranked_res.diff <- res.diff[gene_ranking,]
  counts=as.data.frame(counts(dds))
  
  
  ## distance between curves : used in paper to compare the sense of the upper/lower bias among species
  
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)
  dat0=data.frame(times=times)
   ldc=sapply(row.names(counts),function(gene){
    dat[,levelsfac[1]] <- as.numeric(2 ** (sc.sp1[gene]+sc.fac1[gene,] %*% SDM))
    dat[,levelsfac[2]] <- as.numeric(2 ** (sc.sp2[gene]+sc.fac2[gene,] %*% SDM))
    dat0[,levelsfac[1]] <- as.numeric(2 ** (sc.red.sp1[gene]+sc.red.t[gene,] %*% SDM))
    dat0[,levelsfac[2]]<- as.numeric(2 ** (sc.red.sp2[gene]+sc.red.t[gene,] %*% SDM))
    # distance between curves 
    dc=sum(c((dat[,levelsfac[1]]-dat0[,levelsfac[1]])^2/(dat[,levelsfac[1]]+dat0[,levelsfac[1]]),(dat[,levelsfac[2]]-dat0[,levelsfac[2]])^2/(dat[,levelsfac[2]]+dat0[,levelsfac[2]])))
        return(dc)})
 
   ldcnorm=sapply(row.names(counts),function(gene){
     dat[,levelsfac[1]] <- as.numeric(2 ** (sc.sp1[gene]+sc.fac1[gene,] %*% SDM))
     dat[,levelsfac[2]] <- as.numeric(2 ** (sc.sp2[gene]+sc.fac2[gene,] %*% SDM))
     dat0[,levelsfac[1]] <- as.numeric(2 ** (sc.red.sp1[gene]+sc.red.t[gene,] %*% SDM))
     dat0[,levelsfac[2]]<- as.numeric(2 ** (sc.red.sp2[gene]+sc.red.t[gene,] %*% SDM))
     # distance between curves 
     dc=sum(c((dat[,levelsfac[1]]-dat[,levelsfac[2]])/(dat[,levelsfac[1]]+dat[,levelsfac[2]])))
     return(dc)})
   
  res.diff$dist=ldc
  res.diff$distnormed=ldcnorm
  test_table <- data.frame(res.diff)
  write.table(res.diff, file=paste0(selec,"res.diff.tsv"), sep="\t", quote=F, row.names=T)
  
  return(list(dds,res.diff,sc.sp1,sc.sp2,sc.fac1,sc.fac2,sc.red.sp1,sc.red.sp2,sc.red.t))
}


## to draw profiles 
profile <- function(gene="Pou3f3",facselec="species",dds=dds,selec="md",res.diff,sc.sp1,sc.sp2,sc.fac1,sc.fac2,sc.red.sp1,sc.red.sp2,sc.red.t) {
  df1=colData(dds)
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  k <- counts(dds)[gene,] / sizeFactors(dds)
  SDM <- t(spline_design_matrix(times, intercept=F))
  dat=data.frame(times=times)

  if(facselec=="species")
  {
  dat[,levelsfac[1]] <- as.numeric(2 ** (sc.fac1[gene,] %*% SDM))
  dat[,levelsfac[2]] <- as.numeric(2 ** (sc.fac2[gene,] %*% SDM))
  dat0=data.frame(times=times)
  dat0[,levelsfac[1]] <- as.numeric(2 ** (sc.red.t[gene,] %*% SDM))
  dat0[,levelsfac[2]]<- as.numeric(2 ** (sc.red.t[gene,] %*% SDM))
  }
  if(facselec=="jaw")
  {
    dat[,levelsfac[1]] <- as.numeric(2 ** (sc.sp1[gene]+sc.fac1[gene,] %*% SDM))
    dat[,levelsfac[2]] <- as.numeric(2 ** (sc.sp2[gene]+sc.fac2[gene,] %*% SDM))
    dat0=data.frame(times=times)
    dat0[,levelsfac[1]] <- as.numeric(2 ** (sc.red.sp1[gene]+sc.red.t[gene,] %*% SDM))
    dat0[,levelsfac[2]]<- as.numeric(2 ** (sc.red.sp2[gene]+sc.red.t[gene,] %*% SDM))
  }
  
  dat_melt=melt(dat,id.vars="times")
  dat_melt$type="2 species"
  dat_melt0=melt(dat0,id.vars="times")
  dat_melt0$type="1 curve"
  
  dat0=data.frame(times=df1$timerel,value=k,variable=df1[,fac],rep=df1[,"rep"])
  
  sp_names <- c(`ham` = "hamster",
                `mus` = "mouse")
  to_names <- c(`md` = "lower",
                `mx` = "upper")
  
  models=data.frame(rbind(dat_melt0,dat_melt))
  dat0$jaw=df1[,facselec]
  p=ggplot(dat0, aes(x=times, y=value, group=variable)) + 
    facet_grid(variable~jaw , labeller = labeller(variable=sp_names,jaw=to_names))+theme_bw() +
    geom_point(size=1.5,alpha = 0.2) + theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time (rel)",y="Expr. level (Base Mean)")+
    geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "2 species"="blue"))
  
  save_plot(paste0(id,selec,"rel.png"), ncol=1,nrow=2,
            p ,limitsize = FALSE )   
  return(p)      
}

res_mus_vs_tooth=compute_timeeffect("species","mus","jaw")
res_ham_vs_tooth=compute_timeeffect("species","ham","jaw")

res_md_vs_sp=compute_timeeffect("jaw","md","species")
res_mx_vs_sp=compute_timeeffect("jaw","mx","species")


￼￼## to draw profiles tooth 
profileTooth <- function(gene="Pou3f3",resdds=res_md_vs_sp,facselec="jaw",selec="md") {

      dds=resdds[[1]]
   res.diff=resdds[[2]]
    tsc.fac1=resdds[[3]]
     tsc.fac2=resdds[[4]]
   sc.fac1=resdds[[5]]
   sc.fac2=resdds[[6]]
   sc.red.fac1=resdds[[7]]
   sc.red.fac2=resdds[[8]]
   t.red=resdds[[9]]

  spline_design_matrix <- function(time, intercept = F) {
    tmin <- min(time)
    tmax <- max(time)
    X <- bs(time,
            Boundary.knots = c(tmin,tmax),
            intercept = intercept,
            degree = 3)
    colnames(X) <- paste("bs_", colnames(X), sep = "")
    X
  }
  
  df1=colData(dds)
  df1$timerel[df1$timerel<0]=0
  times = seq(min(df1$timerel),max(df1$timerel),length.out=100)
  k <- counts(dds)[gene,] / sizeFactors(dds)
  SDM <- t(spline_design_matrix(times, intercept=F))
 
  
  if(facselec=="species")
  {
    dat=data.frame(times=times)
    dat[,"md"] <- as.numeric(2 ** (tsc.fac1[gene]+sc.fac1[gene,] %*% SDM))
    dat[,"mx"] <- as.numeric(2 ** (tsc.fac2[gene]+sc.fac2[gene,] %*% SDM))
    dat0=data.frame(times=times)
    dat0[,"md"] <- as.numeric(2 ** (sc.red.fac1[gene]+t.red[gene, ] %*% SDM))
    dat0[,"mx"] <- as.numeric(2 ** (sc.red.fac2[gene]+t.red[gene, ] %*% SDM))
 
  dat_melt=melt(dat,id.vars="times")
  dat_melt0=melt(dat0,id.vars="times")
  dat_melt$type="2 teeth"
  dat_melt0$type="1 curve"
  models=data.frame(rbind(dat_melt0,dat_melt))
  models$sp=selec
  
  dat0=data.frame(times=df1$timerel,value=k,variable=df1[,"jaw"],rep=df1[,"rep"],sp=df1[,"species"])
  
  
  sp_names <- c(ham = "hamster",
                mus = "mouse")
  to_names <- c(md = "lower",
                mx = "upper")
  
  p=ggplot(dat0, aes(x=times, y=value, group=variable)) + 
    facet_grid(sp~variable, labeller=labeller(variable=to_names,sp=sp_names))+theme_bw() +
    geom_point(size=1.5,alpha = 0.2) + theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time",y="Expr. level")+
    geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "2 teeth"="blue"))
  
  
  }
  if(facselec=="jaw")
  {
    # predictions of the 2 curves model
    dat=data.frame(times=times)
    dat[,"ham"] <- as.numeric(2 ** (tsc.fac1[gene]+sc.fac1[gene,] %*% SDM))
    dat[,"mus"] <- as.numeric(2 ** (tsc.fac2[gene]+sc.fac2[gene,] %*% SDM))
    # predictions of the simple model
    dat0=data.frame(times=times)
    dat0[,"ham"] <- as.numeric(2 ** (sc.red.fac1[gene]+t.red[gene, ] %*% SDM))
    dat0[,"mus"] <- as.numeric(2 ** (sc.red.fac2[gene]+t.red[gene, ] %*% SDM))
    
    dat_melt=melt(dat,id.vars="times")
    dat_melt0=melt(dat0,id.vars="times")
    dat_melt$type="2 species"
    dat_melt0$type="1 curve"
    models=data.frame(rbind(dat_melt0,dat_melt))
    
    dat0=data.frame(times=df1$timerel,value=k,variable=df1[,"species"],rep=df1[,"rep"])
    dat0$jaw=selec
    

  
  sp_names <- c(ham = "hamster",
                mus = "mouse")
  to_names <- c(md = "lower",
                mx = "upper")
  
  p=ggplot(dat0, aes(x=times, y=value, group=variable)) + 
    facet_grid(variable~jaw, labeller=labeller(jaw=to_names,variable=sp_names))+theme_bw() +
    geom_point(size=1.5,alpha = 0.2) + theme_bw(base_size = 12)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Dev. time",y="Expr. level")+
    geom_line(data=models,aes(x=times, y=value,col=type,group = type))+
    scale_colour_manual(name="Model",
                        values=c("1 curve"="black",
                                 "2 species"="blue"))
  }
  
  # guides(alpha=FALSE,color=FALSE)
  # ggtitle(paste0(gene," p.value="))
  
  save_plot(paste0(gene,selec,"rel.png"), ncol=1,nrow=2,
            p  )   
  return(p)      
}



plot_mod_sp_vs_Tooth=function(gene,selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
{
 # gene="Pou3f3"
  
  b1=profileTooth(gene,resdds=resdds1,facselec=facsel,selec=selec1)
  b2=profileTooth(gene,resdds=resdds2,facselec=facsel,selec=selec2)
  
  if(facsel=="species")
  {
  prow1 <- plot_grid(
    b2 + theme(legend.position="none"),
    b1 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 2
  )
  legend1 <- get_legend(
    # create some space to the left of the legend
    b1 + theme(legend.box.margin = margin(0, 0, 0, 0))
  )
  
  prow2 <- plot_grid(
    b2 + theme(legend.position="none") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                               axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(strip.background = element_blank(),strip.text.y = element_blank()),
    b1 + theme(legend.position="none") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                               axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(strip.background = element_blank(),strip.text.y = element_blank(),strip.text.x = element_blank()),
    align = 'vh',
    vjust = -1,
    nrow = 2
  )
  
  
  }
  if(facsel!="species")
  {
    prow1 <- plot_grid(
      b1 + theme(legend.position="none"),
      b2 + theme(legend.position="none"),
      align = 'vh',
      labels = c("A", "B"),
      hjust = -1,
      nrow = 1
    )
    legend1 <- get_legend(
      # create some space to the left of the legend
      b1 + theme(legend.box.margin = margin(0, 0, 0, 0))
    )
    
    prow2 <- plot_grid(
      b1 + theme(legend.position="none") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                                 axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
        theme(strip.background = element_blank(),strip.text.y = element_blank()),
      b2 + theme(legend.position="none") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                                 axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
        theme(strip.background = element_blank(),strip.text.y = element_blank()),
      align = 'vh',
      hjust = -1,
      nrow = 1
    )
  }
  
  fig4 <- plot_grid(prow1, legend1, rel_widths = c(1, .22),align = 'bt')
  ggsave(paste0(gene,"_",facsel,"mod2teeth.pdf"),fig4,width = 15, height = 10, units = "cm") 
 
  ggsave(paste0(gene,"_",facsel,"mod2teethb.pdf"),prow2,width = 15, height = 10, units = "cm") 
  
  return(fig4)
}



plot_mod_sp_vs_Tooth("Bmp4",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Bmp4",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")

plot_mod_sp_vs_Tooth("Bmper",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Bmper",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")

plot_mod_sp_vs_Tooth("Fgf3",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Fgf3",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")


plot_mod_sp_vs_Tooth("Fgf10",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Fgf10",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")



plot_mod_sp_vs_Tooth("Bmp6",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Bmp6",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")


plot_mod_sp_vs_Tooth("Bmp2",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Cyp26a1",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")



plot_mod_sp_vs_Tooth("Wif1",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Wif1",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")


plot_mod_sp_vs_Tooth("Dkk1",selec1="mus",selec2="ham",resdds1=res_mus_vs_tooth, resdds2=res_ham_vs_tooth,facsel="species")
plot_mod_sp_vs_Tooth("Dkk1",selec1="md",selec2="mx",resdds1=res_md_vs_sp, resdds2=res_mx_vs_sp,facsel="jaw")


