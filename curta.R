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

Abundance_Index<-function(B,nrep,years,CV,q){
  IA<-matrix(rep(0,length(years)*nrep),nrow=nrep,ncol=length(years))
  for(irep in 1:nrep){
    Index_Abundance<-rep(0,length(years))
    for(i in 1:length(years)){
      v<-(CV*(q)*B[i])^2
      m<-(q)*(B[i])
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

B<-SPM_PT(B_0=B_0,K=K,C=C,p=p,r=r,years=years)
B_SPM<-B

F_SPM <- -log(1-( C/ (( B_SPM+c(B_SPM[-1],B_SPM[length(B_SPM)]))/2)))
F_SPM2 <- C/B_SPM

Bly1<-SPM_PT(B_0=B_SPM[length(B_SPM)],K=K,C=C[length(C)],p=p,r=r,
             years=c(years[length(years)],years[length(years)]+1))
Bly1<-Bly1[length(Bly1)]

B.Bmsy<-Bly1/Bmsy
F.Fmsy<-F_SPM[length(F_SPM)]/Fmsy

cut=30
shape<-B_SPM[1+cut]/Bmsy

n<-p+1

q<-0.000001

SPM_values<-cbind(MSY,Fmsy,Bmsy,K,phi,r,n,shape,B.Bmsy,F.Fmsy,q)
SPM_values


########################## Abundance Index ###########################

nrep=1000
CV=0.5

IA_SPM<-Abundance_Index(B=B,nrep=nrep,years=years,CV=CV,q=q)


############################### Plots ################################

plot(years,B_SPM,type="l",lwd=2,col="red",
     ylab="Biomass",xlab="Time",main="SPM Biomass trajectory")

plot(years,C,type="l",lwd=2,col="red",
     ylab="Captures",xlab="Time",main="Captures trajectory")

IA_SPM_0<-Abundance_Index(B=B_SPM,nrep=1,years=years,CV=0,q=q)
plot(years,IA_SPM_0,type="l",lwd=2,col="red",
     ylim=c((1-CV)*min(IA_SPM_0),(1+CV)*max(IA_SPM_0)),
     ylab="Abundance index",xlab="Time",main="Abundance Index")
for(i in 1:5){
  index<-IA_SPM[i,]
  lines(years,index,lwd=2,col="lightgrey")
}
legend(0.75*length(years)+years[1],0.40*max(IA_SPM_0),c("CV = 0",paste0("CV = ",CV)), 
       lty=c(1,1),col=c("red","grey"))

oldyears<-years
years<-oldyears[-seq(1,cut,by=1)]


#######################################################################
################################ ASPIC ################################
#######################################################################

########################## ASPIC Adjustment ###########################

setwd()

aspic_results<-list()
aspic_noconv<-0
for(irep in 1:nrep){
  inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
            obsI=(IA_SPM[irep,-seq(1,cut,by=1)]), timeI=years,
            ini=list(q=q),
            aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
  waspic(inp,file=paste0(irep,".a7inp"))
  faspic(paste0(irep,".a7inp"),36)
  tt<-file_test("-f",paste0(irep,".rdat"))
  while (tt==FALSE){
    set.seed(round(runif(1,1,1000)))
    aspic_noconv<-aspic_noconv+1/2
    file.remove(paste0(irep,".a7inp"))
    unlink(paste0(irep,".fit"),recursive=T)
    IA_new<-Abundance_Index(B=B_SPM,nrep=1,years=oldyears,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
              obsI=(IA_new[-seq(1,cut,by=1)]), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),36)
    tt<-file_test("-f",paste0(irep,".rdat"))
  }
  res<-raspic(paste0(irep,".fit"))
  errorcode<-as.numeric(res$errorcode)
  while (errorcode!=0) {
    aspic_noconv<-aspic_noconv+1
    file.remove(paste0(irep,".a7inp"))
    file.remove(paste0(irep,".fit"))
    file.remove(paste0(irep,".rdat"))
    IA_new<-Abundance_Index(B=B,nrep=1,years=oldyears,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
              obsI=(IA_new[-seq(1,cut,by=1)]), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),36)
    tt<-file_test("-f",paste0(irep,".rdat"))
    while (tt==FALSE){
      set.seed(round(runif(1,1,1000)))
      aspic_noconv<-aspic_noconv+1/2
      file.remove(paste0(irep,".a7inp"))
      unlink(paste0(irep,".fit"),recursive=T)
      IA_new<-Abundance_Index(B=B_SPM,nrep=1,oldyears=years,CV=CV,q=q)
      inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
                obsI=(IA_new[-seq(1,cut,by=1)]), timeI=years,
                ini=list(q=q),
                aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
      waspic(inp,file=paste0(irep,".a7inp"))
      faspic(paste0(irep,".a7inp"),36)
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
    IA_new<-Abundance_Index(B=B,nrep=1,years=oldyears,CV=CV,q=q)
    inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
              obsI=(IA_new[-seq(1,cut,by=1)]), timeI=years,
              ini=list(q=q),
              aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
    waspic(inp,file=paste0(irep,".a7inp"))
    faspic(paste0(irep,".a7inp"),36)
    tt<-file_test("-f",paste0(irep,".rdat"))
    while (tt==FALSE){
      set.seed(round(runif(1,1,1000)))
      aspic_noconv<-aspic_noconv+1/2
      file.remove(paste0(irep,".a7inp"))
      unlink(paste0(irep,".fit"),recursive=T)
      IA_new<-Abundance_Index(B=B_SPM,nrep=1,years=oldyears,CV=CV,q=q)
      inp<-list(timeC=years, obsC=C[-seq(1,cut,by=1)], 
                obsI=(IA_new[-seq(1,cut,by=1)]), timeI=years,
                ini=list(q=q),
                aspic=list(shape="GENFIT",conditioning="YLD",objfn="SSE"))
      waspic(inp,file=paste0(irep,".a7inp"))
      faspic(paste0(irep,".a7inp"),36)
      tt<-file_test("-f",paste0(irep,".rdat"))
    }
    res<-raspic(paste0(irep,".fit"))
    param_r<-res$estimates$r
  }
  aspic_results[[irep]]<-res
}

aspic_noconv_prop<-aspic_noconv/(nrep+aspic_noconv)
aspic_noconv_prop

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


##################### ASPIC Comparison and Error #####################

ASPIC_estimations<-apply(aspic_estimations,MARGIN=2,mean)
st.dev<-apply(aspic_estimations,MARGIN=2,sd)
aspic_table<-rbind(SPM_values,ASPIC_estimations,st.dev)
rownames(aspic_table)<-c("SPM","ASPIC (media de 100 probas)","ASPIC (st dev)")
colnames(aspic_table)<-c("MSY","Fmsy","Bmsy","K","phi","r","n","shape",
                         "Bly+1/Bmsy","Fly/Fmsy","q")
aspic_table


aspic_RMSE<-matrix(rep(0,nrep*4),ncol=4,nrow=nrep)
colnames(aspic_RMSE)<-c("B","F","B/Bmsy","F/Fmsy")

for(i in 1:nrep){
  predicted<-aspic_states$B0est[i,]
  real<-B_SPM[-seq(1,cut,by=1)]
  aspic_RMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-aspic_states$Fest[i,]
  real<-F_SPM[-seq(1,cut,by=1)]
  aspic_RMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-aspic_states$BBmsy[i,]
  real<-B_SPM[-seq(1,cut,by=1)]/Bmsy
  aspic_RMSE[i,3]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-aspic_states$FFmsy[i,]
  real<- F_SPM[-seq(1,cut,by=1)]/Fmsy
  aspic_RMSE[i,4]<-sqrt(sum((predicted-real)^2)/length(years))  
}

mean_RMSE<-apply(aspic_RMSE,MARGIN=2,mean)
median_RMSE<-apply(aspic_RMSE,MARGIN=2,median)
lq_RMSE<-apply(aspic_RMSE,MARGIN=2,quantile,probs=0.05)
uq_RMSE<-apply(aspic_RMSE,MARGIN=2,quantile,probs=0.95)
aspic_RMSE<-rbind(mean_RMSE,median_RMSE,lq_RMSE,uq_RMSE)
rownames(aspic_RMSE)<-c("media","mediana","lq","uq")
aspic_RMSE


aspic_estimations_RMSE<-matrix(rep(0,ncol(aspic_estimations)),nrow=1)
for(i in 1:ncol(aspic_estimations)){
  real<-aspic_table[1,i]
  aspic_estimations_RMSE[i]<-sqrt(sum((aspic_estimations[,i]-real)^2)/nrep) 
}
colnames(aspic_estimations_RMSE)<-colnames(aspic_estimations)
rownames(aspic_estimations_RMSE)<-"media"
aspic_estimations_RMSE


aspic_MAPE<-matrix(rep(0,nrep*4),ncol=4,nrow=nrep)
colnames(aspic_MAPE)<-c("B","F","B/Bmsy","F/Fmsy")

for(i in 1:nrep){
  predicted<-aspic_states$B0est[i,]
  real<-B_SPM[-seq(1,cut,by=1)]
  aspic_MAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$Fest[i,]
  real<-F_SPM[-seq(1,cut,by=1)]
  aspic_MAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$BBmsy[i,]
  real<-B_SPM[-seq(1,cut,by=1)]/Bmsy
  aspic_MAPE[i,3]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$FFmsy[i,]
  real<- F_SPM[-seq(1,cut,by=1)]/Fmsy
  aspic_MAPE[i,4]<-1/length(years)*(sum(abs((predicted-real)/real))) 
}

mean_MAPE<-apply(aspic_MAPE,MARGIN=2,mean)
median_MAPE<-apply(aspic_MAPE,MARGIN=2,median)
lq_MAPE<-apply(aspic_MAPE,MARGIN=2,quantile,probs=0.05)
uq_MAPE<-apply(aspic_MAPE,MARGIN=2,quantile,probs=0.95)
aspic_MAPE<-rbind(mean_MAPE,median_MAPE,lq_MAPE,uq_MAPE)
rownames(aspic_MAPE)<-c("media","mediana","lq","uq")
aspic_MAPE

aspic_estimations_MAPE<-matrix(rep(0,ncol(aspic_estimations)),nrow=1)
for(i in 1:ncol(aspic_estimations)){
  real<-aspic_table[1,i]
  aspic_estimations_MAPE[i]<- 1/nrep*(sum(abs((real-aspic_estimations[,i])/real)))
}
colnames(aspic_estimations_MAPE)<-colnames(aspic_estimations)
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
for(i in 1:nrep){

  inp1=list(timeC=years, obsC=C[-seq(1,cut,by=1)],
            obsI=IA_SPM[i,-seq(1,cut,by=1)], timeI=years)
  res<- fit.spict(inp1)
  convergence<-res$opt$convergence

  while (convergence!=0) {
    spict_noconv<-spict_noconv+1
    IA_new<-Abundance_Index(B=B,nrep=1,years=oldyears,CV=CV,q=q)
    inp1=list(timeC=years, obsC=C[-seq(1,cut,by=1)],
              obsI=IA_new[-seq(1,cut,by=1)], timeI=years)
    res<- fit.spict(inp1)
    convergence<-res$opt$convergence
  }

  inp <- res$inp

  iB <- get.par("logB", res, exp = TRUE)[seq(1,length(inp$time),by=16),2]
  B_spict<-rbind(B_spict,iB)
  iF <- get.par("logFs", res, exp = TRUE)[seq(1,length(inp$time),by=16), 2]
  F_spict<-rbind(F_spict,iF)

  MSY<-get.par("logMSY",res,exp=TRUE)[,2]
  Fmsy<-get.par("logFmsy",res,exp=TRUE)[,2]
  Bmsy<-get.par("logBmsy",res,exp=TRUE)[,2]
  K<-get.par("logK",res,exp=TRUE)[,2]
  phi<-Bmsy/K
  r<-get.par("logr",res,exp=TRUE)[,2]
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

  iestimations<-cbind(MSY,Fmsy,Bmsy,K,phi,r,n,shape,B.Bmsy,F.Fmsy,q)
  spict_estimations<-rbind(spict_estimations,iestimations)
  cat(100*i/nrep,"%\n")

}

spict_noconv_prop<-spict_noconv/(nrep+spict_noconv)
spict_noconv_prop

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


##################### SPiCT Comparison and Error #####################

SPICT_estimations<-apply(spict_estimations,MARGIN=2,mean)
st.dev<-apply(spict_estimations,MARGIN=2,sd)
spict_table<-rbind(SPM_values,SPICT_estimations,st.dev)
rownames(spict_table)<-c("SPM","SPiCT (media de 100 probas)","SPiCT (st dev)")
colnames(spict_table)<-c("MSY","Fmsy","Bmsy","K","phi","r","n","shape",
                         "Bly+1/Bmsy","Fly/Fmsy","q")
spict_table


spict_RMSE<-matrix(rep(0,nrep*4),ncol=4,nrow=nrep)
colnames(spict_RMSE)<-c("B","F","B/Bmsy","F/Fmsy")

for(i in 1:nrep){
  predicted<-spict_states$B[i,]
  real<-B_SPM[-seq(1,cut,by=1)]
  spict_RMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-spict_states$F[i,]
  real<-F_SPM[-seq(1,cut,by=1)]
  spict_RMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-spict_states$BR[i,]
  real<-B_SPM[-seq(1,cut,by=1)]/Bmsy.r
  spict_RMSE[i,3]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-spict_states$FR[i,]
  real<- F_SPM[-seq(1,cut,by=1)]/Fmsy
  spict_RMSE[i,4]<-sqrt(sum((predicted-real)^2)/length(years))  
}

mean_RMSE<-apply(spict_RMSE,MARGIN=2,mean)
median_RMSE<-apply(spict_RMSE,MARGIN=2,median)
lq_RMSE<-apply(spict_RMSE,MARGIN=2,quantile,probs=0.05)
uq_RMSE<-apply(spict_RMSE,MARGIN=2,quantile,probs=0.95)
spict_RMSE<-rbind(mean_RMSE,median_RMSE,lq_RMSE,uq_RMSE)
rownames(spict_RMSE)<-c("media","mediana","lq","uq")
spict_RMSE


spict_estimations_RMSE<-matrix(rep(0,ncol(spict_estimations)),nrow=1)
for(i in 1:ncol(spict_estimations)){
  real<-spict_table[1,i]
  spict_estimations_RMSE[i]<-sqrt(sum((spict_estimations[,i]-real)^2)/nrep) 
}
colnames(spict_estimations_RMSE)<-colnames(spict_estimations)
rownames(spict_estimations_RMSE)<-"media"
spict_estimations_RMSE


spict_MAPE<-matrix(rep(0,nrep*4),ncol=4,nrow=nrep)
colnames(spict_MAPE)<-c("B","F","B/Bmsy","F/Fmsy")

for(i in 1:nrep){
  predicted<-spict_states$B[i,]
  real<-B_SPM[-seq(1,cut,by=1)]
  spict_MAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$F[i,]
  real<-F_SPM[-seq(1,cut,by=1)]
  spict_MAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$BR[i,]
  real<-B_SPM[-seq(1,cut,by=1)]/Bmsy.r
  spict_MAPE[i,3]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$FR[i,]
  real<- F_SPM[-seq(1,cut,by=1)]/Fmsy
  spict_MAPE[i,4]<-1/length(years)*(sum(abs((predicted-real)/real))) 
}

mean_MAPE<-apply(spict_MAPE,MARGIN=2,mean)
median_MAPE<-apply(spict_MAPE,MARGIN=2,median)
lq_MAPE<-apply(spict_MAPE,MARGIN=2,quantile,probs=0.05)
uq_MAPE<-apply(spict_MAPE,MARGIN=2,quantile,probs=0.95)
spict_MAPE<-rbind(mean_MAPE,median_MAPE,lq_MAPE,uq_MAPE)
rownames(spict_MAPE)<-c("media","mediana","lq","uq")
spict_MAPE

spict_estimations_MAPE<-matrix(rep(0,ncol(spict_estimations)),nrow=1)
for(i in 1:ncol(spict_estimations)){
  real<-spict_table[1,i]
  spict_estimations_MAPE[i]<- 1/nrep*(sum(abs((real-spict_estimations[,i])/real)))
}
colnames(spict_estimations_MAPE)<-colnames(spict_estimations)
rownames(spict_estimations_MAPE)<-"media"
spict_estimations_MAPE


######################################################################
########################## ASPIC vs SPiCT ############################
######################################################################

com_table<-rbind(SPM_values,ASPIC_estimations,SPICT_estimations)
rownames(com_table)<-c("SPM","ASPIC","SPiCT")
colnames(com_table)<-c("MSY","Fmsy","Bmsy","K","phi","r","n","shape",
                       "Bly+1/Bmsy","Fly/Fmsy","q")
taboa_comparativa<-com_table
taboa_comparativa
