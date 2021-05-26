library(expm)
library(maxLik)
library(ggplot2)
library(reshape2)
library(plyr)
library(cowplot)
library("readxl")
library(MASS)
library(mgcv)
library(ggeffects)
#remotes::install_github("gavinsimpson/gratia")
library("gratia")
library("mgcv")



### MODELLING THE CUSP PATTERNING IN MOUSE AND HAMSTER


###   MOUSE MOLARS

# Reading the raw data. 
cuspmxM <- read_excel("cusp_fgf4_mouse.xlsx", sheet = "mx2020")
cuspmdM <- read_excel("cusp_fgf4_mouse.xlsx", sheet = "md2020")
cuspmdM=cuspmdM[,c("age","averageweight","totalcuspids","stage")]
names(cuspmdM)=c("age","weight","total","stage")
cuspmxM=cuspmxM[,c("age","averageweight","totalcuspids","stage")]
names(cuspmxM)=c("age","weight","total","stage")


###   HAMSTER MOLARS
# Reading the raw data. 
cuspmxH <- read_excel("cusp_fgf4_hamster.xlsx", sheet = "mx2020")
cuspmdH <- read_excel("cusp_fgf4_hamster.xlsx", sheet = "md2020")
cuspmxH=cuspmxH[,c("age","averageweight","totalcuspids","stage")]
names(cuspmxH)=c("age","weight","total","stage")
cuspmdH=cuspmdH[,c("age","averageweight","totalcuspids","stage")]
names(cuspmdH)=c("age","weight","total","stage")



### Sample for weight/time relationships
df_weight_time=read.table("weightTime_hamstermouse.txt",h=T)

### Samples for RNAseq phenotype
metadata=read.table("metadata-samplesRNAseq.txt",h=T)




log_weight_by_species = function(data=df_weight_time,species="mus",confidence="yes"){
  df = data[data$species == species,]
  if(confidence=="yes")
  {
    df = df[df$confidence=="yes",]
  }
  df$log_weight= log(df$weight)
  stage=df$stage
  log_weight=df$log_weight
  lm_out = lm(stage~log_weight)
  plot(lm_out$fitted, lm_out$resid)
  plot(df$log_weight,(df$stage),pch=20,xlab="log_weight",ylab=paste("time"),main=paste(species,"time/ log(weight)"))
  abline(lm_out)
  #print(summary(lm_out))
  return(lm_out)
}


boxcox_by_species = function(data, species,confidence="yes"){
  df = data[data$species == species,]
  if(confidence=="yes")
  {
    df = df[df$confidence=="yes",]
  }
  bc = boxcox(df$weight~df$stage, lambda=seq(-1,1,.001))
  r=data.frame(cbind(bc$x,bc$y))
  lambda = r[order(-r$X2),][1,1]
  weight_lambda = df$weight^lambda
  stage = df$stage
  lm_out = lm(stage~weight_lambda)
  plot(lm_out$fitted, lm_out$resid)
  plot(df$weight^lambda,(df$stage),pch=20,xlab=paste("weight^",lambda),ylab="time",main=paste(species,"time/ BC weight"))
  abline(lm_out)
  #print(summary(lm_out))
  return(list(lambda,lm_out))
}


gam_by_species = function(data, species,confidence="yes"){
  
  df = data[data$species == species,]
  if(confidence=="yes")
  {
    df = df[df$confidence=="yes",]
  }
  bc = boxcox(df$weight~df$stage, lambda=seq(-1,1,.001))
  r=data.frame(cbind(bc$x,bc$y))
  lambda = r[order(-r$X2),][1,1]
  weight_lambda = df$weight^lambda
  stage = df$stage
  gam_out = gam(stage~ s(weight_lambda, bs="cr"))
  lm_out = gam(stage~ weight_lambda)
  
  anova(lm_out, gam_out, test="Chisq")
  pred.gam <- ggpredict(gam_out, terms = c("weight_lambda"))  # this gives overall predictions for the model
  # Plot the predictions 
 gam=(ggplot(pred.gam) + 
   geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                     fill = "lightgrey", alpha = 0.5) +  # error band
         geom_point(data = df,                      # adding the raw data (scaled values)
                    aes(x = weight_lambda, y = stage)) + 
         labs(x = "Weight (lambda)", y = "stage", 
              title = "Weight and stage model") + 
         theme_minimal()
     )
 ggsave(file=paste0("gam_",species,"_",confidence,".pdf"),gam)
 
 pred.gam$weight=(pred.gam$x)^(1/lambda)
 gam=(ggplot(pred.gam) + 
        geom_line(aes(x = weight, y = predicted)) +          # slope
        geom_ribbon(aes(x = weight, ymin = predicted - std.error, ymax = predicted + std.error), 
                    fill = "lightgrey", alpha = 0.5) +  # error band
        geom_point(data = df,                      # adding the predicted data (scaled values)
                   aes(x = weight, y = stage),color="black",alpha=0.1) + 
        labs(x = "Weight", y = "stage", 
             title = "Weight and stage model") + 
        theme_minimal()
 )
 ggsave(file=paste0("gam_rev_",species,"_",confidence,".pdf"),gam)
 
   return(list(lambda,gam_out))
}







predict_species_withRNAseq = function(species="mus", data=df_weight_time, newdataP=cuspmdM,data_name="cuspmdM",confidence="yes",metadatatab=metadata){
  log_model = log_weight_by_species(data, species,confidence=confidence)
  bc_out = boxcox_by_species(data, species,confidence=confidence)
  lambda = bc_out[[1]]
  bc_model = bc_out[[2]]
  gam_out = gam_by_species(data, species,confidence=confidence)
  gam_model = gam_out[[2]]
  
  newdataP=newdataP[!is.na(newdataP$weight),]
  newdat_log <- data.frame(
    log_weight = log(newdataP$weight))
  newdat_bc = data.frame(
    weight_lambda = newdataP$weight^(lambda)) 
  
  x = predict(log_model,newdata = newdat_log,se.fit=TRUE, interval = "prediction", level = 0.95)
  y = predict(bc_model,newdata = newdat_bc,se.fit=TRUE, interval = "prediction", level = 0.95)
  z = confint(gam_model, parm = "s(weight_lambda)", type = "confidence", newdata=newdat_bc)
  
  
  tooth=ifelse(grepl("md",data_name),"md","mx")
  newdataRNASeq=metadatatab[metadatatab$espece==species&metadatatab$machoire==tooth,c("echantillon","stade","replicat","weight")]
  
  rnaseqdat_log <- data.frame(
  log_weight = log(newdataRNASeq$weight))
  rnaseqdat_bc = data.frame(
  weight_lambda = newdataRNASeq$weight^(lambda)) 
  
  xr = predict(log_model,newdata = rnaseqdat_log,se.fit=TRUE, interval = "prediction", level = 0.95)
  yr = predict(bc_model,newdata = rnaseqdat_bc,se.fit=TRUE, interval = "prediction", level = 0.95)
  zr = confint(gam_model, parm = "s(weight_lambda)", type = "confidence", newdata=rnaseqdat_bc)

  
  plot(x$fit[,"fit"],y$fit[,"fit"],xlab="predictions log model",ylab="predictions bc model")
  points(xr$fit[,"fit"],yr$fit[,"fit"],pch=20,col="red")
  abline(0,1)
  
  plot(x$fit[,"fit"],(z$est+gam_model$coefficients["(Intercept)"]),xlab="predictions log model",ylab="predictions gam model")
  points(xr$fit[,"fit"],(zr$est+gam_model$coefficients["(Intercept)"]),pch=20,col="red")
  abline(0,1)
  
  
  zr$fit=zr$est+gam_model$coefficients["(Intercept)"]
  zr$lower1=zr$lower+gam_model$coefficients["(Intercept)"]
  zr$upper1=zr$upper+gam_model$coefficients["(Intercept)"]
  
  z$fit=z$est+gam_model$coefficients["(Intercept)"]
  z$lower1=z$lower+gam_model$coefficients["(Intercept)"]
  z$upper1=z$upper+gam_model$coefficients["(Intercept)"]
  
  
  pred.gam.cusp = z
  pred.gam.rnaseq = zr
  

  pred.gam.cusp$weight=newdataP$weight
  pred.gam.rnaseq$weight=newdataRNASeq$weight
  pred.gam.cusp$stage=newdataP$age
  pred.gam.rnaseq$stage=newdataRNASeq$stade
  
  
  gam=(ggplot(pred.gam.cusp) + 
         geom_line(aes(x = weight, y = fit)) +          # slope
         geom_ribbon(aes(x = weight, ymin = lower1, ymax = upper1), 
                     fill = "lightgrey", alpha = 0.5) +  # error band
         
         geom_point(data = pred.gam.cusp,                      # adding the predicted data (scaled values)
                    aes(x = weight, y = fit),color="dark grey",alpha=0.8) + 
         
         geom_point(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                    aes(x = weight, y = fit),color="dark red",alpha=0.8) + 
         geom_errorbar(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                       aes(x = weight, ymin = lower1, ymax = upper1), colour="black", width=.1) +
         geom_point(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                    aes(x = weight, y = stage),color="dark orange",alpha=0.8)+
         labs(x = "Weight", y = "stage", 
              title = "Weight and stage model") + 
         theme_minimal()
  )
  ggsave(file=paste0("gam_rev_rnaseq_",data_name,"_",confidence,".pdf"),gam)
  
  pdf(file=paste0("predictions_RNAseq",data_name,"_",confidence,".pdf"))
  plot(pred.gam.rnaseq$weight,pred.gam.rnaseq$fit,xlab="weight",ylab="predictions",pch=20,col="dark red")
  points(pred.gam.rnaseq$weight,yr$fit[,"fit"],pch=20,col="blue")
  points(pred.gam.rnaseq$weight,pred.gam.rnaseq$stage,pch=20,col="dark orange")
  points(pred.gam.rnaseq$weight,xr$fit[,"fit"],pch=20,col="dark green")
  legend(x="bottomright",col=c("dark red","blue","dark green","dark orange"),legend=c("boxCox+GAM","boxCox","log","dpc"),pch=20)
  dev.off()
  
  yr=data.frame(yr$fit)
  names(yr)=paste0(names(yr),"_boxCox")
  xr=data.frame(xr$fit)
  names(xr)=paste0(names(xr),"_log")
  pred.gam.rnaseq=zr[,c("fit","lower1","upper1")]
  names(pred.gam.rnaseq)=paste0(names(pred.gam.rnaseq),"_GAM")
  
  newdataRNASeq=data.frame(cbind(newdataRNASeq,yr,xr,pred.gam.rnaseq))
  write.table(file=paste0("predictions_rnaseq_",data_name,"_",confidence,".txt"),newdataRNASeq)
  
  y=data.frame(y$fit)
  names(y)=paste0(names(y),"_boxCox")
  x=data.frame(x$fit)
  names(x)=paste0(names(x),"_log")
  z=z[,c("fit","lower1","upper1")]
  names(z)=paste0(names(z),"_GAM")
  
  
  newdataP=data.frame(cbind(newdataP,x,y,z))
  newdataP$total2=NULL
  newdataP$total2=newdataP$total
  newdataP$total2[newdataP$total%in%c("PEK/SEK","no","SEK faible","PEK/SEK faible","SEK","PEK/SEK")]=1
  newdataP$total2[newdataP$total=="BUD"|newdataP$total=="bud"]=-1
  newdataP$total2[newdataP$total=="PEK"]=0
  write.table(file=paste0("predictions_cusp_",data_name,"_",confidence,".txt"),newdataP)
  
  ### to predict for all possible weights
  
  newdataTweight=seq(from=min(na.omit(newdataP$weight)),to=max(na.omit(newdataP$weight)),by=1)
  newdatT_log <- data.frame(
    log_weight = log(newdataTweight))
  newdatT_bc = data.frame(
    weight_lambda = newdataTweight^(lambda)) 
  xt = predict(log_model,newdata = newdatT_log,se.fit=TRUE, interval = "prediction", level = 0.95)
  yt = predict(bc_model,newdata = newdatT_bc,se.fit=TRUE, interval = "prediction", level = 0.95)
 zt = confint(gam_model, parm = "s(weight_lambda)", type = "confidence", newdata=newdatT_bc)
  
  zt$fit=zt$est+gam_model$coefficients["(Intercept)"]
  zt$lower1=zt$lower+gam_model$coefficients["(Intercept)"]
  zt$upper1=zt$upper+gam_model$coefficients["(Intercept)"]
   
  yt=data.frame(yt$fit)
  names(yt)=paste0(names(yt),"_boxCox")
  xt=data.frame(xt$fit)
  names(xt)=paste0(names(xt),"_log")
  zt=zt[,c("fit","lower1","upper1")]
  names(zt)=paste0(names(zt),"_GAM")
  
  newdataT=data.frame(cbind(newdataTweight,xt,yt,zt))
  write.table(file=paste0("predictions_timeline_",data_name,"_",confidence,".txt"),newdataT)
  
  return(newdataP)
}



cuspmdM_all=predict_species_withRNAseq(species="mus", data=df_weight_time, newdataP=cuspmdM, data_name="cuspmdM",confidence="all" )
cuspmxM_all=predict_species_withRNAseq(species="mus", data=df_weight_time, newdataP=cuspmxM, data_name="cuspmxM",confidence="all" )
cuspmxH_all=predict_species_withRNAseq(species="ham",data=df_weight_time, newdataP=cuspmxH, data_name="cuspmxH",confidence="all" )
cuspmdH_all=predict_species_withRNAseq(species="ham",data=df_weight_time, newdataP=cuspmdH, data_name="cuspmdH",confidence="all" )


cuspmxM_s=cuspmxM_all[,c("fit_GAM","total2")]
cuspmdM_s=cuspmdM_all[,c("fit_GAM","total2")]
cuspmxH_s=cuspmxH_all[,c("fit_GAM","total2")]
cuspmdH_s=cuspmdH_all[,c("fit_GAM","total2")]
cuspmdH_s$total2[cuspmdH_s$total2==7]=6


cuspmxM_b=cuspmxM_all[,c("fit_boxCox","total2")]
cuspmdM_b=cuspmdM_all[,c("fit_boxCox","total2")]
cuspmxH_b=cuspmxH_all[,c("fit_boxCox","total2")]
cuspmdH_b=cuspmdH_all[,c("fit_boxCox","total2")]
cuspmdH_b$total2[cuspmdH_b$total2==7]=6



cuspmxM_s$total2=as.numeric(cuspmxM_s$total2)
cuspmxM_s=na.omit(cuspmxM_s)
cuspmdM_s$total2=as.numeric(cuspmdM_s$total2)
cuspmdM_s=na.omit(cuspmdM_s)
cuspmdH_s$total2=as.numeric(cuspmdH_s$total2)
cuspmdH_s=na.omit(cuspmdH_s)
cuspmxH_s$total2=as.numeric(cuspmxH_s$total2)
cuspmxH_s=na.omit(cuspmxH_s)
names(cuspmxM_s)=c("time","cusps")
names(cuspmdM_s)=c("time","cusps")
names(cuspmxH_s)=c("time","cusps")
names(cuspmdH_s)=c("time","cusps")


cuspmxM_b$total2=as.numeric(cuspmxM_b$total2)
cuspmxM_b=na.omit(cuspmxM_b)
cuspmdM_b$total2=as.numeric(cuspmdM_b$total2)
cuspmdM_b=na.omit(cuspmdM_b)
cuspmdH_b$total2=as.numeric(cuspmdH_b$total2)
cuspmdH_b=na.omit(cuspmdH_b)
cuspmxH_b$total2=as.numeric(cuspmxH_b$total2)
cuspmxH_b=na.omit(cuspmxH_b)
names(cuspmxM_b)=c("time","cusps")
names(cuspmdM_b)=c("time","cusps")
names(cuspmxH_b)=c("time","cusps")
names(cuspmdH_b)=c("time","cusps")




### models for cusp addition


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

# rep :  k of the gamma law
# returns the log likelihood

calculLikTx=function(v,init,pattern,rep,datacusp){
  rep2=c(rep(1:length(rep),rep),length(rep)+1)
  Q=fQparam(v,pattern,rep)
  probs=apply(datacusp,1,function(x){
    t=as.numeric(x[1])-init #bidouille sinon pb trop de time à 0 cusp
    mat=expm(t*Q)
    cusp=as.numeric(x[2])+2 #pour transformer -1 en 1, 0 en 2 etc
    return(sum(mat[1,rep2==cusp]))})
  lik=sum(log(probs))
  return(lik)
}

# for plotting with bars centered on max probability 

plotmodel2=function(rep=c(5,5,5,5,5,5,5,5,5),pattern=c(1:9),data=cuspmxM_s,model=tauxMx9Tx,title='ppmx',end=19,init=13.5)
{
  rep2=c(rep(1:length(rep),rep),length(rep)+1)
  t=seq(init,end,by=0.1)
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

  ggsave(paste(title,".png",sep=""),pmmx,width = 20, height = 15, units = "cm") 
  ggsave(paste(title,".pdf",sep=""),pmmx,width = 20, height = 15, units = "cm") 
  
  save(model,file=paste(title,".RData",sep=""))
  
  return(pmmx)
}


### MODELS FOR UPPER MOLAR
# Model with 9 different rates
# maxLik maximizes the log likelihood (the biggest lnL the best)
# init = first stage in dpc (13.5 for mouse)

init=min(cuspmxM_s$time)-0.00001
mini=min(cuspmxM_s$time)
maxi=max(cuspmxM_s$time)
pattern9=c(1:9)
rep=c(5,5,5,5,5,5,5,5,5)
MusMx9Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=pattern9,rep=rep,datacusp=cuspmxM_s))},
  start=1+abs(rnorm(length(rep))),
  method="BFGS",
  iterlim = 100,
  constraints=list(ineqA=rbind(diag(1,length(unique(pattern9))),diag(-1,length(unique(pattern9)))),
                   ineqB=c(rep(0,length(unique(pattern9))),rep(40,length(unique(pattern9))))))
m9=plotmodel2(rep=c(5,5,5,5,5,5,5,5,5),pattern=pattern9,model=MusMx9Tx,title='pmmx9center',end=maxi,init=mini)



# Model with 4 different rates

rep=c(5,5,5,5,5,5,5,5,5)
pattern4=c(4,1,3,2,1,1,2,2,1)
MusMx4Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=pattern4,rep=rep,datacusp=cuspmxM_s))},
                 start=1+abs(rnorm(length(unique(pattern4)))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(pattern4))),diag(-1,length(unique(pattern4)))),
                                 ineqB=c(rep(0,length(unique(pattern4))),rep(40,length(unique(pattern4))))))
m4=plotmodel2(rep=c(5,5,5,5,5,5,5,5,5),pattern=pattern4,model=MusMx4Tx,title='pmmx4center',end=maxi,init=mini)


# Model with 3 different rates
rep=c(5,5,5,5,5,5,5,5,5)
pattern3=c(2,1,3,2,1,1,2,2,1)
MusMx3Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=pattern3,rep=rep,datacusp=cuspmxM_s))},
                 start=1+abs(rnorm(length(unique(pattern3)))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(pattern3))),diag(-1,length(unique(pattern3)))),
                                 ineqB=c(rep(0,length(unique(pattern3))),rep(40,length(unique(pattern3))))))
m3=plotmodel2(rep=c(5,5,5,5,5,5,5,5,5),pattern=pattern3,model=MusMx3Tx,title='pmmx3center',end=maxi,init=mini)

# Model with 2 different rates
pattern2=c(2,1,1,2,1,1,2,2,1)
MusMx2Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=pattern2,rep=rep,datacusp=cuspmxM_s))},
                 start=1+abs(rnorm(length(unique(pattern2)))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(pattern2))),diag(-1,length(unique(pattern2)))),
                                 ineqB=c(rep(0,length(unique(pattern2))),rep(40,length(unique(pattern2))))))
m2=plotmodel2(rep=c(5,5,5,5,5,5,5,5,5),pattern=pattern2,model=MusMx2Tx,title='pmmx2center',end=maxi,init=mini)

# Model with 1 rate
pattern1=c(1,1,1,1,1,1,1,1,1)
MusMx1Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=pattern1,rep=rep,datacusp=cuspmxM_s))},
                 start=1+abs(rnorm(length(unique(pattern1)))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(pattern1))),diag(-1,length(unique(pattern1)))),
                                 ineqB=c(rep(0,length(unique(pattern1))),rep(40,length(unique(pattern1))))))
m1=plotmodel2(rep=c(5,5,5,5,5,5,5,5,5),pattern=pattern1,model=MusMx1Tx,title='pmmx1center',end=maxi,init=mini)


proba9Tx4Tx= 1-pchisq(2*(MusMx9Tx$maximum-MusMx4Tx$maximum),df=5)
proba9Tx3Tx= 1-pchisq(2*(MusMx9Tx$maximum-MusMx3Tx$maximum),df=6)
proba4Tx3Tx= 1-pchisq(2*(MusMx4Tx$maximum-MusMx3Tx$maximum),df=1)
proba3Tx2Tx= 1-pchisq(2*(MusMx3Tx$maximum-MusMx2Tx$maximum),df=1)
proba2Tx1Tx= 1-pchisq(2*(MusMx2Tx$maximum-MusMx1Tx$maximum),df=1)





### MODELS FOR LOWER MOLAR IN MOUSE

init=min(cuspmdM_s$time)-0.00001
mini=min(cuspmdM_s$time)
maxi=max(cuspmdM_s$time)
rep=c(5,5,5,5,5,5,5)

# Model with 8 different rates
patternmd8=c(1:7)
MusMd8Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=patternmd8,rep=rep,datacusp=cuspmdM_s))},
                 start=1+abs(rnorm(length(rep))),
                 method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternmd8))),diag(-1,length(unique(patternmd8)))),
                                 ineqB=c(rep(0,length(unique(patternmd8))),rep(40,length(unique(patternmd8))))))

# Model with 3 different rates
patternmd3=c(2,2,1,1,3,3,3)
MusMd3Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=patternmd3,rep=rep,datacusp=cuspmdM_s))},
                 start=1+abs(rnorm(length(unique(patternmd3)))),
                method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternmd3))),diag(-1,length(unique(patternmd3)))),
                ineqB=c(rep(0,length(unique(patternmd3))),rep(40,length(unique(patternmd3))))))

  

# Model with 2 different rates
patternmd2=c(2,2,1,1,2,2,2)
MusMd2Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=patternmd2,rep=rep,datacusp=cuspmdM_s))},
                 start=1+abs(rnorm(length(unique(patternmd2)))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternmd2))),diag(-1,length(unique(patternmd2)))),
                                 ineqB=c(rep(0,length(unique(patternmd2))),rep(40,length(unique(patternmd2))))))

# Model with 1 rate
init=min(cuspmdM_s$time)
mini=min(cuspmdM_s$time)
maxi=max(cuspmdM_s$time)
rep=c(5,5,5,5,5,5,5)
patternmd1=c(1,1,1,1,1,1,1)
MusMd1Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=patternmd1,rep=rep,datacusp=cuspmdM_s))},
                 start=1+abs(rnorm(length(patternmd1))),
                 method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternmd1))),diag(-1,length(unique(patternmd1)))),
                                 ineqB=c(rep(0,length(unique(patternmd1))),rep(40,length(unique(patternmd1))))))



md8=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternmd8,model=MusMd8Tx,title='pmmd8center',data=cuspmdM_s,end=maxi,init=mini)
md2=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternmd2,model=MusMd2Tx,title='pmmd2center',data=cuspmdM_s,end=maxi,init=mini)
md3=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternmd3,model=MusMd3Tx,title='pmmd3center',data=cuspmdM_s,end=maxi,init=mini)


proba8Tx3Tx= 1-pchisq(2*(MusMd8Tx$maximum-MusMd3Tx$maximum),df=5)
proba8Tx2Tx= 1-pchisq(2*(MusMd8Tx$maximum-MusMd2Tx$maximum),df=6)
proba3Tx2Tx= 1-pchisq(2*(MusMd3Tx$maximum-MusMd2Tx$maximum),df=1)



###### MOLARS IN HAMSTER

##### LOWER HAM
init=min(cuspmdH_s$time)-0.00001
mini=min(cuspmdH_s$time)
mini=12
init=mini
maxi=max(cuspmdH_s$time)
rep=c(5,5,5,5,5,5,5)
patternhd1=c(1,1,1,1,1,1,1)
HamMd1Tx=maxLik(function(v){return(calculLikTx(v,init=init,pattern=patternhd1,rep=rep,datacusp=cuspmdH_s))},
                start=1+abs(rnorm(length(unique(patternhd1)))),
                method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhd1))),diag(-1,length(unique(patternhd1)))),ineqB=c(rep(0,length(unique(patternhd1))),rep(40,length(unique(patternhd1))))))
mdh1=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhd1,model=HamMd1Tx,title='phmd1center',data=cuspmdH_s,end=maxi,init=mini)


patternhd2=c(2,1,1,1,2,2,2)
HamMd2Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhd2,rep=rep,datacusp=cuspmdH_s))},
                start=1+abs(rnorm(length(unique(patternhd2)))),
                method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhd2))),diag(-1,length(unique(patternhd2)))),ineqB=c(rep(0,length(unique(patternhd2))),rep(40,length(unique(patternhd2))))))
mdh2=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhd2,model=HamMd2Tx,title='phmd2center1',data=cuspmdH_s,end=maxi,init=mini)


patternhd3=c(3,2,2,2,1,1,1)
HamMd3Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhd3,rep=rep,datacusp=cuspmdH_s))},
                start=1+abs(rnorm(length(unique(patternhd3)))),
                method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhd3))),diag(-1,length(unique(patternhd3)))),ineqB=c(rep(0,length(unique(patternhd3))),rep(40,length(unique(patternhd3))))))
mdh3=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhd3,model=HamMd3Tx,title='phmd3center',data=cuspmdH_s,end=maxi,init=mini)


patternhd7=c(1,2,3,4,5,6,6)
HamMd7Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhd7,rep=rep,datacusp=cuspmdH_s))},
                start=1+abs(rnorm(length(unique(patternhd7)))),
                method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhd7))),
                                             diag(-1,length(unique(patternhd7)))),
                                 ineqB=c(rep(0,length(unique(patternhd7))),
                                         rep(40,length(unique(patternhd7))))))
mdh7=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhd7,model=HamMd7Tx,title='phmd7center',data=cuspmdH_s,end=maxi,init=mini)




proba7Tx3Tx= 1-pchisq(2*(HamMd7Tx$maximum-HamMd3Tx$maximum),df=4)
proba7Tx2Tx= 1-pchisq(2*(HamMd7Tx$maximum-HamMd2Tx$maximum),df=5)
proba2Tx1Tx= 1-pchisq(2*(HamMd2Tx$maximum-HamMd1Tx$maximum),df=1)
proba3Tx2Tx= 1-pchisq(2*(HamMd3Tx$maximum-HamMd2Tx$maximum),df=1)

proba7Tx3Tx
proba7Tx2Tx
proba2Tx1Tx
proba3Tx2Tx
proba7Tx3Tx
proba7Tx3Tx


### ## UPPER HAM

mini=min(cuspmxH_s$time)
maxi=max(cuspmxH_s$time)
mini=12
init=mini

rep=c(5,5,5,5,5,5,5)
patternhx1=c(1,1,1,1,1,1,1)
HamMx1Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhx1,rep=rep,datacusp=cuspmxH_s))},
                start=1+abs(rnorm(length(unique(patternhx1)))),
                method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhx1))),
                                             diag(-1,length(unique(patternhx1)))),
                                 ineqB=c(rep(0,length(unique(patternhx1))),
                                         rep(40,length(unique(patternhx1))))))



patternhx2=c(2,1,1,2,2,2,2)
HamMx2Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhx2,rep=rep,datacusp=cuspmxH_s))},
                start=1+abs(rnorm(length(unique(patternhx2)))),
                method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhx2))),diag(-1,length(unique(patternhx2)))),ineqB=c(rep(0,length(unique(patternhx2))),rep(40,length(unique(patternhx2))))))
mxh2=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhx2,model=HamMx2Tx,title='phmx2center',data=cuspmxH_s,end=maxi,init=mini)



patternhx3=c(2,1,1,3,3,3,3)
HamMx3Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhx3,rep=rep,datacusp=cuspmxH_s))},
                start=1+abs(rnorm(length(unique(patternhx3)))),
                method="BFGS",
                iterlim = 40,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhx3))),diag(-1,length(unique(patternhx3)))),ineqB=c(rep(0,length(unique(patternhx3))),rep(40,length(unique(patternhx3))))))
mxh3=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhx3,model=HamMx3Tx,title='phmx3center',data=cuspmxH_s,end=maxi,init=mini)




patternhx7=c(1,2,3,4,5,6,7)
HamMx7Tx=maxLik(function(v){return(calculLikTx(v,init=11.8,pattern=patternhx7,rep=rep,datacusp=cuspmxH_s))},
                start=1+abs(rnorm(length(unique(patternhx7)))),
                method="BFGS",
                iterlim = 200,
                constraints=list(ineqA=rbind(diag(1,length(unique(patternhx7))),
                                             diag(-1,length(unique(patternhx7)))),
                                 ineqB=c(rep(0,length(unique(patternhx7))),
                                         rep(40,length(unique(patternhx7))))))
mxh7=plotmodel2(rep=c(5,5,5,5,5,5,5),pattern=patternhx7,model=HamMx7Tx,title='phmx7center',data=cuspmxH_s,end=maxi,init=mini)



proba7Tx3Tx= 1-pchisq(2*(HamMx7Tx$maximum-HamMx3Tx$maximum),df=4)
proba2Tx1Tx= 1-pchisq(2*(HamMx2Tx$maximum-HamMx1Tx$maximum),df=1)
proba3Tx2Tx= 1-pchisq(2*(HamMx3Tx$maximum-HamMx2Tx$maximum),df=1)





#### FIGURE DRAWING WITH BEST MODELS


models=plot_grid(m3, mxh3, md3, mdh3, labels = c('A', 'B', 'C', 'D'), label_size = 16, align = "v")
modelsannotate=plot_grid(m3, mxh3, md3, mdh3, labels = c('Mus Mx', 'Ham Mx', 'Mus Md', 'Ham Md'), label_size = 16, hjust=-1)
ggsave("modelsannotate.png",modelsannotate,width = 30, height = 30, units = "cm") 
ggsave("models.png",models,width = 30, height = 30, units = "cm") 


ggsave("modelsannotate.pdf",modelsannotate,width = 30, height = 30, units = "cm") 
ggsave("models.pdf",models,width = 30, height = 30, units = "cm") 


models=plot_grid(m9, mxh7, md8, mdh7, labels = c('A', 'B', 'C', 'D'), label_size = 16, align = "v")
modelsannotate=plot_grid(m9, mxh7, md8, mdh7,labels = c('Mus Mx', 'Ham Mx', 'Mus Md', 'Ham Md'), label_size = 16, hjust=-1)
ggsave("modelsannotatecomplexe.png",modelsannotate,width = 30, height = 30, units = "cm") 
ggsave("modelscomplexe.png",models,width = 30, height = 30, units = "cm") 

ggsave("modelsannotatecomplexe.pdf",modelsannotate,width = 30, height = 30, units = "cm") 
ggsave("modelscomplexe.pdf",models,width = 30, height = 30, units = "cm") 

