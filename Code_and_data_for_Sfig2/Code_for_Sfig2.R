library(cowplot)
library(expm)
library(plyr)
library(ggplot2)


# for plotting with bars centered on max probability 


fQparam=function(v,p,rep)
{
  Q=matrix(data=0,nrow=sum(rep)+1,ncol=sum(rep)+1)
  r2=rep(1:length(rep),rep)
  z=sapply(1:length(r2),function(x){
    Q[x,x]<<-(-v[p[r2[x]]])
    Q[x,x+1]<<-v[p[r2[x]]]
  })
  return(Q)
}

plotmodel_bars=function(rep=c(5,5,5,5,5,5,5,5,5),pattern=c(1:9),data=cuspmxM_s,model=tauxMx9Tx,title='Mouse upper',end=19,init=13.5)
{
  rep2=c(rep(1:length(rep),rep),length(rep)+1)
  t=seq(init,end,by=0.01)
  tm=t-init
  # mp contains the probability to be in each state at each time point
  mp=sapply(tm,function(tp){
    Q=fQparam(model$estimate,pattern,rep)
    mat=expm(tp*Q)
    p=tapply(mat[1,],rep2,sum)
    return(p)
  })
  tmax=apply(mp,1,which.max)
  x=rep(model$estimate[pattern],rep)
  cc1=c(0,cumsum(1/x))
  rep2=c(rep(1:length(rep),rep),length(rep)+1)
  le=sapply(1:max(rep2),function(i){ma=ifelse(i<max(rep2),min(cc1[rep2==i+1]),max(cc1));
  duree=ma-min(cc1[rep2==i])})
  mini=sapply(1:max(rep2),function(i){t[tmax][i]-le[i]/2})
  maxi=sapply(1:max(rep2),function(i){t[tmax][i]+le[i]/2})
  
  df=data.frame(time=data$time,cusp=data$cusps)
  df1=data.frame(cusp=-1:max(data$cusps),
                 modelstart=mini,
                 modelend=maxi,
                 cuspend=factor(c(0:max(data$cusps),max(data$cusps))), 
                 tmax=t[tmax])
  df1$cusp2=factor(df1$cusp)
  if(min(df$cusp)==-1)
  {
    df1$cusp2=revalue(df1$cusp2, c("-1"="Bud","0"="PEK"))
    df1$cusp2=as.character( df1$cusp2)
    df$cusp2=factor(df$cusp)
    df$cusp2=revalue(df$cusp2, c("-1"="Bud","0"="PEK"))
  }
  if(min(df$cusp)==0)
  {
    df1$cusp2=revalue(df1$cusp2, c("0"="PEK"))
    df1$cusp2=as.character( df1$cusp2)
    df$cusp2=factor(df$cusp)
    df$cusp2=revalue(df$cusp2, c("0"="PEK"))
  }
  
  
  pmmx=ggplot(df, aes(x=time, y=cusp2)) + 
    geom_point() + 
    scale_y_discrete("Cusps",limits=c("Bud","PEK",1:8))+
    scale_x_continuous("Dev. time") + #, limits=c(init,end))+
    geom_segment(aes( x = modelstart, y = cusp2, xend = modelend, yend = cusp2), data = df1,color="blue",size=5,alpha=0.2)

#  ggsave(paste(title,".png",sep=""),pmmx,width = 20, height = 15, units = "cm") 
#  ggsave(paste(title,".pdf",sep=""),pmmx,width = 20, height = 15, units = "cm") 
  
  return(pmmx)
}


# HAMSTER UPPER
model_hamupper<-readRDS("datamodel_HamUpper_for_figS2.rds")
cusp_hamupper<-read.table("datacusp_HamUpper_for_figS2.txt",sep="\t",h=T)
mini=min(cusp_hamupper$time)-0.00001
maxi=max(cusp_hamupper$time)
init=mini
pattern5=c(1:5)
rep=c(5,5,5,5,5)
plot_hamupper=plotmodel_bars(rep=rep,pattern=pattern5,model=model_hamupper,title='Ham Upper',data=cusp_hamupper,end=maxi,init=mini)


# HAMSTER LOWER
model_hamlower<-readRDS("datamodel_HamLower_for_figS2.rds")
cusp_hamlower<-read.table("datacusp_HamLower_for_figS2.txt",sep="\t",h=T)
mini=min(cusp_hamlower$time)-0.00001
maxi=max(cusp_hamlower$time)
init=mini
plot_hamlower=plotmodel_bars(rep=rep,pattern=pattern5,model=model_hamlower,title='Ham Lower',data=cusp_hamlower,end=maxi,init=mini)


# MOUSE LOWER
model_muslower<-readRDS("datamodel_MusLower_for_figS2.rds")
cusp_muslower<-read.table("datacusp_MusLower_for_figS2.txt",sep="\t",h=T)
init=min(cusp_muslower$time)
mini=min(cusp_muslower$time)
maxi=max(cusp_muslower$time)
plot_muslower=plotmodel_bars(rep=rep,pattern=pattern5,model=model_muslower,title='Mus Lower',data=cusp_muslower,end=maxi,init=mini)

# MOUSE UPPER
model_musupper<-readRDS("datamodel_MusUpper_for_figS2.rds")
cusp_musupper<-read.table("datacusp_MusUpper_for_figS2.txt",sep="\t",h=T)
init=min(cusp_musupper$time)
mini=min(cusp_musupper$time)
maxi=max(cusp_musupper$time)
plot_musupper=plotmodel_bars(rep=rep,pattern=pattern5,model=model_musupper,title='Mus Upper',data=cusp_musupper,end=maxi,init=mini)



# ASSEMBLY

# this corresponds to the development age of samples with lower and upper time boundaries
metadataRNAseq<-read.table("metadata_RNAseq.txt",h=T)
MusZero=mean(metadataRNAseq[metadataRNAseq$echantillon%in%c("mus_AXR_MAOS","mus_BGF_MAOS"),]$fit_GAM)
HamZero=mean(metadataRNAseq[metadataRNAseq$echantillon%in%c("ham_BCR_AOSW"),]$fit_GAM)
MusTen=mean(metadataRNAseq[metadataRNAseq$echantillon%in%c("mus_AXR_MQOS","mus_BGF_MQOS"),]$fit_GAM)
HamTen=mean(metadataRNAseq[metadataRNAseq$echantillon%in%c("ham_BCR_QOSW","ham_BGE_QOSW"),]$fit_GAM)


modelsannotate=plot_grid(plot_hamupper+ theme_classic() + geom_vline( xintercept=HamZero)+ geom_vline( xintercept=HamTen),
                         plot_hamlower+ theme_classic()+ geom_vline( xintercept=HamZero)+ geom_vline( xintercept=HamTen),
                         plot_musupper+ theme_classic()+ xlim(13.5,18.25) + geom_vline( xintercept=MusZero)+ geom_vline( xintercept=MusTen), 
                         plot_muslower+ theme_classic()+ xlim(13.5,18.21)+ geom_vline( xintercept=MusZero)+ geom_vline( xintercept=MusTen),
                         labels = c('Ham Mx', 'Ham Md', 'Mus Mx', 'Mus Md'), label_size = 16, hjust=-1)
ggsave("figS2.pdf",modelsannotate,width = 30, height = 30, units = "cm") 


