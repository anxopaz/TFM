library(Rfishpop)
library(connectASPIC)
library(lamW)
library(spict)


years=seq(2020,2100,1)
ctrPop<-list(years=years,niter=1,N0=15000,ages=0:15,minFage=2,
             maxFage=5,tc=0.5,seed=123)
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
Mvec=c(1,0.6,0.5,0.4,0.35,0.35,0.3,rep(0.3,9))
M<-matrix(rep(Mvec,number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages
ctrBio<-list(M=M,CV_M=0, L_inf=20, t0=0, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), 
             b=3,a50_Mat=-4, ad_Mat=-0.2,CV_Mat=0)
ctrSEL<-list(type="cte", par=list(cte=1),CV_SEL=0)
f=matrix(rep(0.5,number_years),ncol=number_years,nrow=1,byrow=TRUE)
ctrFish<-list(f=f,ctrSEL=ctrSEL)
a_BH=15000; b_BH=50; CV_REC_BH=0
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

RE<-BYR.eq(Pop.Mod,0,3,3,c(FALSE,1),Method="mean",par=NULL)
N_eq<-RE$N
rf=RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1,plot=FALSE)
fmsy=rf$F_msy[,1,1];fmsy

ctrPop<-list(years=years,niter=1,N0=N_eq,ages=0:15,minFage=2,
             maxFage=5,ts=0,tc=0.5,seed=123)

f_1=c(0.001,rep(fmsy,15),
      seq(fmsy,6*fmsy,length.out=15),
      seq(6*fmsy,0.01,length.out=5),
      seq(0.01,fmsy,length.out=4),
      rep(fmsy,number_years-40))

f=matrix(f_1,ncol=number_years,nrow=1,byrow=TRUE)
ctrFish<-list(f=f,ctrSEL=ctrSEL)
CV_REC_BH=0
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

C.r=Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,1]

rf=RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1,plot=FALSE)
MSY.r<-rf$F_msy[,6,1]
Fmsy.r<-rf$F_msy[,2,1]
Bmsy.r<-rf$F_msy[,7,1]
eq0<-BYR.eq(Pop.Mod,0,3,3,F,Method="mean",par=NULL)
K.r<-eq0$DPM[,7,1]


############################# Functions ##############################

SPM_PT<-function(B_0,K,C,p,r,years){
  B<-rep(0,length(years))
  B[1]<-B_0
  for(i in 1:(length(years)-1)) {
    B[i+1]<-B[i]+(r/p)*B[i]*(1-(B[i]/K)^p)-C[i]
  }
  return(B)
}

SPM_PT2<-function(B_0,K,C,p,r,years,cv){
  B<-rep(0,length(years))
  B[1]<-B_0
  for(i in 1:(length(years)-1)) {
    m<-r
    v<-(cv*m)^2
    mu<-log(m^2/sqrt(m^2+v))
    sigma<-sqrt(log(v/m^2+1))
    r_aux=rlnorm(1,mu,sigma)
    B[i+1]<-B[i]+(r_aux/p)*B[i]*(1-(B[i]/K)^p)-C[i]
  }
  return(B)
}

Abundance_Index<-function(B,nrep,years,CV,q){
  IA<-matrix(rep(0,length(years)*nrep),nrow=nrep,ncol=length(years))
  for(irep in 1:nrep){
    Index_Abundance<-rep(0,length(years))
    for(i in 1:length(years)){
      v<-(CV*(q)*B[irep,i])^2
      m<-(q)*(B[irep,i])
      mu<-log(m^2/sqrt(m^2+v))
      sigma<-sqrt(log(v/m^2+1))
      Index_Abundance[i]<-stats::rlnorm(1, meanlog =mu, sdlog =sigma)
    }
    IA[irep,]<-Index_Abundance
  }
  return(IA)
}


###################### Surplus Production Model ######################

MSY<-MSY.r
Bmsy<-Bmsy.r
Fmsy<-Fmsy.r
B_0<-K.r
K<-K.r

phi=Bmsy/K.r
p=(lambertWm1(phi*log(phi))-log(phi))/log(phi)
r=(Fmsy*p)/(1-phi^p)

C<-C.r
cv_p<-0.3

nrep=1000

B<-matrix(rep(0,length(years)*nrep),nrow=nrep)
for(i in 1:nrep){
  B2<-SPM_PT2(B_0,K,C,p,r,years,cv_p)
  while (is.na(B2[length(years)])==TRUE){
    B2<-SPM_PT2(B_0,K,C,p,r,years,cv_p)
  }
  while (B2[length(years)]<=0){
    B2<-SPM_PT2(B_0,K,C,p,r,years,cv_p)
  }
  B[i,]<-B2
}
B_SPM<-B

B0=SPM_PT(B_0,K,C,p,r,years)

F_SPM <- C/ (( B0+c(B0[-1],B0[length(B0)]))/2)

Bly1<-SPM_PT(B_0=B0[length(B0)],K=K,C=C[length(C)],p=p,r=r,
             years=c(years[length(years)],years[length(years)]+1))
Bly1<-Bly1[length(Bly1)]

B.Bmsy<-Bly1/Bmsy
F.Fmsy<-F_SPM[length(F_SPM)]/Fmsy
shape<-B0[1]/Bmsy
n<-p+1

q<-0.000001

SPM_values<-cbind(MSY,Fmsy,Bmsy,K,phi,r,n,shape,B.Bmsy,F.Fmsy,q)
SPM_values


########################## Abundance Index ###########################

CV=0.3

IA_SPM<-Abundance_Index(B=B,nrep=nrep,years=years,CV=CV,q=q)


############################### Plots ################################
plot(years,B0,type="l",ylab="Biomasa",xlab="Ano",main="Biomasa")
for(i in 1:100){
  lines(years,B[i,],col="lightgrey")}
lines(years,B0,lwd=2,col="red")
legend(x="bottomright",c("sen erro de proceso","con erro"), 
       lty=c(1,1),col=c("red","grey"))

plot(years,C,type="l",lwd=2,col="red",
     ylab="Capturas",xlab="Ano",main="Traxectoria das Capturas")

IA_SPM_0<-B0*q
plot(years,IA_SPM_0,type="l",lwd=2,col="red",
     ylim=c((1-CV)*min(IA_SPM_0),(1+CV)*max(IA_SPM_0)),
     ylab="Índice de biomasa",xlab="Tempo",main="Índice de abundancia")
for(i in 1:15){
  index<-IA_SPM[i,]
  lines(years,index,lwd=2,col="lightgrey")
}
lines(years,IA_SPM_0,type="l",lwd=2,col="red")
legend(0.75*length(years)+years[1],0.40*max(IA_SPM_0),c("CV = 0",paste0("CV = ",CV)), 
       lty=c(1,1),col=c("red","grey"))


#######################################################################
################################ ASPIC ################################
#######################################################################

########################## ASPIC Adjustment ###########################

setwd("~/Desktop/MTE/TFM/ASPIC/")

aspic_results<-list()
aspic_noconv<-0
Ba<-matrix(rep(0,length(years)*nrep),nrow=nrep)

for(irep in 1:nrep){
  inp<-list(timeC=years, obsC=C, 
            obsI=(IA_SPM[irep,]), timeI=years,
            ini=list(q=q),
            aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
  B2<-B[irep,]
  waspic(inp,file=paste0(irep,".a7inp"))
  faspic(paste0(irep,".a7inp"),t=30)
  tt<-file_test("-f",paste0(irep,".rdat"))
  while (tt==FALSE){
    aspic_noconv<-aspic_noconv+1/2
    file.remove(paste0(irep,".a7inp"))
    unlink(paste0(irep,".fit"),recursive=T)
    B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    while (sum(is.na(B2)!=0)){
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    }
    IA_new<-Abundance_Index(B=as.matrix(B2),nrep=1,years=years,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C, 
              obsI=(IA_new), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),t=30)
    tt<-file_test("-f",paste0(irep,".rdat"))
  }
  res<-raspic(paste0(irep,".fit"))
  errorcode<-as.numeric(res$errorcode)
  while (errorcode!=0) {
    aspic_noconv<-aspic_noconv+1
    file.remove(paste0(irep,".a7inp"))
    file.remove(paste0(irep,".fit"))
    file.remove(paste0(irep,".rdat"))
    B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    while (sum(is.na(B2)!=0)){
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    }
    IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C, 
              obsI=(IA_new), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),t=30)
    tt<-file_test("-f",paste0(irep,".rdat"))
    while (tt==FALSE){
      aspic_noconv<-aspic_noconv+1
      file.remove(paste0(irep,".a7inp"))
      unlink(paste0(irep,".fit"),recursive=T)
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
      while (sum(is.na(B2)!=0)){
        B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
      }
      IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
      inp<-list(timeC=years, obsC=C, 
                obsI=(IA_new), timeI=years,
                ini=list(q=q),
                aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
      waspic(inp,file=paste0(irep,".a7inp"))
      faspic(paste0(irep,".a7inp"),t=30)
      tt<-file_test("-f",paste0(irep,".rdat"))
    }
    res<-raspic(paste0(irep,".fit"))
    errorcode<-as.numeric(res$errorcode)
  }
  param_r<-res$estimates$r
  while (is.na(param_r)==TRUE){
    aspic_noconv<-aspic_noconv+1
    file.remove(paste0(irep,".a7inp"))
    file.remove(paste0(irep,".fit"))
    file.remove(paste0(irep,".rdat"))
    B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    while (sum(is.na(B2)!=0)){
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    }
    IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C, 
              obsI=(IA_new), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),t=30)
    tt<-file_test("-f",paste0(irep,".rdat"))
    while (tt==FALSE){
      aspic_noconv<-aspic_noconv+1
      file.remove(paste0(irep,".a7inp"))
      unlink(paste0(irep,".fit"),recursive=T)
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
      while (sum(is.na(B2)!=0)){
        B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
      }
      IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
      inp<-list(timeC=years, obsC=C, 
                obsI=(IA_new), timeI=years,
                ini=list(q=q),
                aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
      waspic(inp,file=paste0(irep,".a7inp"))
      faspic(paste0(irep,".a7inp"),t=30)
      tt<-file_test("-f",paste0(irep,".rdat"))
    }
    
    res<-raspic(paste0(irep,".fit"))
    param_r<-res$estimates$r
  }
  param_n<-res$estimates$n
  while (param_n==1){
    aspic_noconv<-aspic_noconv+1
    file.remove(paste0(irep,".a7inp"))
    file.remove(paste0(irep,".fit"))
    file.remove(paste0(irep,".rdat"))
    B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C, 
              obsI=(IA_new), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),t=30)
    tt<-file_test("-f",paste0(irep,".rdat"))
    while (tt==FALSE){
      aspic_noconv<-aspic_noconv+1
      file.remove(paste0(irep,".a7inp"))
      unlink(paste0(irep,".fit"),recursive=T)
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
      IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
      inp<-list(timeC=years, obsC=C, 
                obsI=(IA_new), timeI=years,
                ini=list(q=q),
                aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
      waspic(inp,file=paste0(irep,".a7inp"))
      faspic(paste0(irep,".a7inp"),t=30)
      tt<-file_test("-f",paste0(irep,".rdat"))
    }
    
    res<-raspic(paste0(irep,".fit"))
    param_n<-res$estimates$n
  }
  Ba[irep,]<-B2
  aspic_results[[irep]]<-res
}

aspic_noconv_prop<-aspic_noconv/(nrep+aspic_noconv)
cat("Probabilidade de non converxencia = ",aspic_noconv_prop)
names(aspic_results) <- seq(1,nrep,by=1)

nstates <- colnames(aspic_results[[1]]$states)
aspic_states <- vector("list", length(nstates))
names(aspic_states) <- nstates

for(i in 1:nrep){
  
  mat <- aspic_results[[i]]$states
  
  for(istates in nstates){
    aspic_states[[istates]] <- rbind(aspic_states[[istates]],mat[,istates])
  }
  
}

for(istates in nstates){
  colnames(aspic_states[[istates]]) <- years
  rownames(aspic_states[[istates]])<-seq(1,nrep,by=1)
}

aspic_estimations=NULL
for(i in 1:nrep){
  iestimations<-cbind((aspic_results[[i]]$estimates))
  aspic_estimations<-rbind(aspic_estimations,iestimations)
}


####################### ASPIC real observations ######################

names<-c("MSY","Bmsy","BT/Bmsy","K","n")

RMS<-NULL
Brms<-NULL
BTBrms<-NULL

for (i in 1:nrep) {
  biomasa_t<-Ba[i,-length(Ba[i,])]
  biomasa_t1<-c(Ba[i,-1])
  producion<-biomasa_t1-biomasa_t+C[-length(C)]
  index<-which.max(producion)
  iRMS<-producion[index]
  iBrms<-biomasa_t[index]
  RMS<-c(RMS,iRMS)
  Brms<-c(Brms,iBrms)
  iBTBrms<-Ba[i,length(Ba[i,])]/iBrms
  BTBrms<-c(BTBrms,iBTBrms)
}

Kre<-c(rep(K,nrep))
nre<-c(rep(n,nrep))
aspic_real<-cbind(RMS,Brms,BTBrms,Kre,nre)
colnames(aspic_real)<-names



##################### ASPIC Comparison and Error #####################

aspic_estimations_plot<-aspic_estimations[,c(1,3,9,4,7)]
colnames(aspic_estimations_plot)<-names

aspic_estimations_mean<-apply(aspic_estimations_plot,MARGIN=2,mean)
aspic_estimations_sd<-apply(aspic_estimations_plot,MARGIN=2,sd)


aspic_RMSE<-matrix(rep(0,nrep*2),ncol=2,nrow=nrep)
colnames(aspic_RMSE)<-c("B","B/Bmsy")

for(i in 1:nrep){
  predicted<-aspic_states$B0est[i,]
  real<-Ba[i,]
  aspic_RMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-aspic_states$BBmsy[i,]
  real<-Ba[i,]/Brms[i]
  aspic_RMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
}

mean_RMSE<-apply(aspic_RMSE,MARGIN=2,mean)
median_RMSE<-apply(aspic_RMSE,MARGIN=2,median)
lq_RMSE<-apply(aspic_RMSE,MARGIN=2,quantile,probs=0.05)
uq_RMSE<-apply(aspic_RMSE,MARGIN=2,quantile,probs=0.95)
aspic_RMSE<-rbind(mean_RMSE,median_RMSE,lq_RMSE,uq_RMSE)
rownames(aspic_RMSE)<-c("media","mediana","lq","uq")
aspic_RMSE

aspic_estimations_RMSE<-matrix(rep(0,ncol(aspic_estimations_plot)),nrow=1)
for(i in 1:ncol(aspic_estimations_plot)){
  real<-aspic_real[,i]
  aspic_estimations_RMSE[i]<-sqrt(sum((aspic_estimations_plot[,i]-real)^2)/nrep) 
}
colnames(aspic_estimations_RMSE)<-names
rownames(aspic_estimations_RMSE)<-"media"
aspic_estimations_RMSE

aspic_MAPE<-matrix(rep(0,nrep*2),ncol=2,nrow=nrep)
colnames(aspic_MAPE)<-c("B","B/Bmsy")

for(i in 1:nrep){
  predicted<-aspic_states$B0est[i,]
  real<-Ba[i,]
  aspic_MAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$BBmsy[i,]
  real<-Ba[i,]/Brms[i]
  aspic_MAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
}

mean_MAPE<-apply(aspic_MAPE,MARGIN=2,mean)
median_MAPE<-apply(aspic_MAPE,MARGIN=2,median)
lq_MAPE<-apply(aspic_MAPE,MARGIN=2,quantile,probs=0.05)
uq_MAPE<-apply(aspic_MAPE,MARGIN=2,quantile,probs=0.95)
aspic_MAPE<-rbind(mean_MAPE,median_MAPE,lq_MAPE,uq_MAPE)
rownames(aspic_MAPE)<-c("media","mediana","lq","uq")
aspic_MAPE

aspic_estimations_MAPE<-matrix(rep(0,ncol(aspic_estimations_plot)),nrow=1)
for(i in 1:ncol(aspic_estimations_plot)){
  real<-aspic_real[,i]
  aspic_estimations_MAPE[i]<- 1/nrep*(sum(abs((real-aspic_estimations_plot[,i])/real)))
}
colnames(aspic_estimations_MAPE)<-names
rownames(aspic_estimations_MAPE)<-"media"
aspic_estimations_MAPE



#######################################################################
################################ SPiCT ################################
#######################################################################

######################### SPiCT Adjustment ############################


B_spict=NULL; F_spict=NULL; BR_spict<-NULL; FR_spict<-NULL
spict_rep_states<-list()
spict_estimations=NULL
spict_noconv<-0
Bs<-matrix(rep(0,length(years)*nrep),nrow=nrep)

for(i in 1:nrep){

  inp1=list(timeC=years, obsC=C, obsI=IA_SPM[i,], timeI=years)
  res<- fit.spict(inp1)
  convergence<-res$opt$convergence
  B2<-B[i,]

  while (convergence!=0) {
    spict_noconv<-spict_noconv+1
    B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    while (sum(is.na(B2)!=0)){
      B2<-matrix(SPM_PT2(B_0,K,C,p,r,years,cv_p),nrow=1)
    }
    IA_new<-Abundance_Index(B=B2,nrep=1,years=years,CV=CV,q=q)
    inp1=list(timeC=years, obsC=C, obsI=IA_new, timeI=years)
    res<- fit.spict(inp1)
    convergence<-res$opt$convergence
  }

  Bs[i,]<-B2
  inp <- res$inp

  iB <- get.par("logB", res, exp = TRUE)[seq(1,length(inp$time),by=16),2]
  B_spict<-rbind(B_spict,iB)
  iF <- get.par("logFs", res, exp = TRUE)[seq(1,length(inp$time),by=16), 2]
  F_spict<-rbind(F_spict,iF)

  MSY<-get.par("logMSY",res,exp=TRUE)[,2]
  Fmsy<-get.par("logFmsy",res,exp=TRUE)[,2]
  Bmsy<-get.par("logBmsy",res,exp=TRUE)[,2]
  K.s<-get.par("logK",res,exp=TRUE)[,2]
  phi<-Bmsy/K
  r.s<-get.par("logr",res,exp=TRUE)[,2]
  n<-get.par("p",res)[,2] + 1
  shape<-iB[1]/Bmsy
  B.Bmsy<-get.par("logBlBmsy",res,exp=TRUE)[,2]
  F.Fmsy<-get.par("logFlFmsy",res,exp=TRUE)[,2]
  q<-get.par("logq",res,exp=TRUE)[,2]

  iBR<-iB/Bmsy
  iFR<-iF/Fmsy
  BR_spict<-rbind(BR_spict,iBR)
  FR_spict<-rbind(FR_spict,iFR)

  spict_rep_states[[i]]<-rbind(iB,iF,iBR,iFR)
  rownames(spict_rep_states[[i]])<-c("B","F","BR","FR")

  iestimations<-cbind(MSY,Fmsy,Bmsy,K=K.s,phi,r.s,n,shape,B.Bmsy,F.Fmsy,q)
  spict_estimations<-rbind(spict_estimations,iestimations)

  cat(100*i/nrep,"%\n")

}

spict_noconv_prop<-spict_noconv/(nrep+spict_noconv)
cat("Probabilidade de non converxencia = ",spict_noconv_prop)

spict_states<-list(B_spict=B_spict,F_spict=F_spict,BR_spict=BR_spict,FR_spict=FR_spict)
spict_names<-c("B","F","BR","FR")
names(spict_states)<-spict_names

for(inames in spict_names){
  rownames(spict_states[[inames]])<-1:nrep
}

spict_states_plus2<-spict_states

spict_states<-list()
for (i in 1:4){ 
  spict_states[[i]]<-spict_states_plus2[[i]][,-c(length(years)+1,
                                                 length(years)+2)]
}


names(spict_states)<-names(spict_states_plus2)
colnames(spict_states$F)<-years
colnames(spict_states$FR)<-years


####################### SPiCT real observations ######################

names<-c("MSY","Bmsy","BT/Bmsy","K","n")

RMS2<-NULL
Brms2<-NULL
BTBrms2<-NULL

for (i in 1:nrep) {
  biomasa_t<-Bs[i,-length(Bs[i,])]
  biomasa_t1<-c(Bs[i,-1])
  producion<-biomasa_t1-biomasa_t+C[-length(C)]
  index<-which.max(producion)
  iRMS<-producion[index]
  iBrms<-biomasa_t[index]
  RMS2<-c(RMS2,iRMS)
  Brms2<-c(Brms2,iBrms)
  iBTBrms<-Ba[i,length(Bs[i,])]/iBrms
  BTBrms2<-c(BTBrms2,iBTBrms)
}

spict_real<-cbind(RMS2,Brms2,BTBrms2,Kre,nre)
colnames(spict_real)<-names



##################### SPiCT Comparison and Error #####################

spict_estimations_plot<-spict_estimations[,c(1,3,9,4,7)]
colnames(spict_estimations_plot)<-names

spict_estimations_mean<-apply(spict_estimations_plot,MARGIN=2,mean)
spict_estimations_sd<-apply(spict_estimations_plot,MARGIN=2,sd)


spict_RMSE<-matrix(rep(0,nrep*2),ncol=2,nrow=nrep)
colnames(spict_RMSE)<-c("B","B/Bmsy")

for(i in 1:nrep){
  predicted<-spict_states$B[i,]
  real<-Bs[i,]
  spict_RMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years))
  predicted<-spict_states$BR[i,]
  real<-Bs[i,]/Brms2[i]
  spict_RMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
}

mean_RMSE<-apply(spict_RMSE,MARGIN=2,mean)
median_RMSE<-apply(spict_RMSE,MARGIN=2,median)
lq_RMSE<-apply(spict_RMSE,MARGIN=2,quantile,probs=0.05)
uq_RMSE<-apply(spict_RMSE,MARGIN=2,quantile,probs=0.95)
spict_RMSE<-rbind(mean_RMSE,median_RMSE,lq_RMSE,uq_RMSE)
rownames(spict_RMSE)<-c("media","mediana","lq","uq")
spict_RMSE

spict_estimations_RMSE<-matrix(rep(0,ncol(spict_estimations_plot)),nrow=1)
for(i in 1:ncol(spict_estimations_plot)){
  real<-spict_real[,i]
  spict_estimations_RMSE[i]<-sqrt(sum((spict_estimations_plot[,i]-real)^2)/nrep) 
}
colnames(spict_estimations_RMSE)<-names
rownames(spict_estimations_RMSE)<-"media"
spict_estimations_RMSE

spict_MAPE<-matrix(rep(0,nrep*2),ncol=2,nrow=nrep)
colnames(spict_MAPE)<-c("B","B/Bmsy")

for(i in 1:nrep){
  predicted<-spict_states$B[i,]
  real<-Bs[i,]
  spict_MAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$BR[i,]
  real<-Bs[i,]/Brms2[i]
  spict_MAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
}

mean_MAPE<-apply(spict_MAPE,MARGIN=2,mean)
median_MAPE<-apply(spict_MAPE,MARGIN=2,median)
lq_MAPE<-apply(spict_MAPE,MARGIN=2,quantile,probs=0.05)
uq_MAPE<-apply(spict_MAPE,MARGIN=2,quantile,probs=0.95)
spict_MAPE<-rbind(mean_MAPE,median_MAPE,lq_MAPE,uq_MAPE)
rownames(spict_MAPE)<-c("media","mediana","lq","uq")
spict_MAPE

spict_estimations_MAPE<-matrix(rep(0,ncol(spict_estimations_plot)),nrow=1)
for(i in 1:ncol(spict_estimations_plot)){
  real<-spict_real[,i]
  spict_estimations_MAPE[i]<- 1/nrep*(sum(abs((real-spict_estimations_plot[,i])/real)))
}
colnames(spict_estimations_MAPE)<-names
rownames(spict_estimations_MAPE)<-"media"
spict_estimations_MAPE

