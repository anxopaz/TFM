library(Rfishpop)
library(lamW)
library(ggplot2)
library(gridExtra)


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

SPM_PT<-function(B_0,K,C,p,r,years){
  B<-rep(0,length(years))
  B[1]<-B_0
  for(i in 1:(length(years)-1)) {
    B[i+1]<-B[i]+(r/p)*B[i]*(1-(B[i]/K)^p)-C[i]
  }
  return(B)
}

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


plot1<-data.frame(Biomasa=B_SPM,Anos=years)
plot2<-data.frame(Capturas=C,Anos=years)
plot3<-data.frame(Esforzo=f_1,Anos=years,fmsy)

p1<-ggplot(plot1, aes(x=Anos, y=Biomasa)) +
  geom_line() +
  labs(x="Anos",y="Biomasa (t)",title="Traxectoria da Biomasa") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=16,face="bold",hjust=0.5))

p2<-ggplot(plot2, aes(x=Anos, y=Capturas)) +
  geom_line() +
  labs(x="Anos",y="Capturas (t)", title="Serie das Capturas") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=16,face="bold",hjust=0.5))

p3<-ggplot(plot3, aes(x=Anos, y=Esforzo)) +
  geom_line() +
  labs(x="Anos",y="Mortalidade por pesca", title="Mortalidade por pesca") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=16,face="bold",hjust=0.5)) + 
  geom_hline(yintercept = fmsy,linetype = "dashed",col="red")


Idades<-seq(0,15,by=1)
LS<-Length_VB(L_inf=20,k=0.3,ages=Idades,t0=0)
plot4<-data.frame(Idades,LS)
p4<-ggplot(plot4, aes(x=Idades, y=LS)) +
  geom_line() +
  geom_point() +
  labs(x="Idades",y="Talla do Stock (cm)",title="Talla do stock") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))


WS<-Weight(L=LS , a=6*10^(-6),b=3)
plot5<-data.frame(LS,WS)
p5<-ggplot(plot5, aes(x=LS, y=WS)) +
  geom_line() +
  geom_point() +
  labs(x="Talla do Stock (cm)",y="Peso do Stock (Kg)",title="Relación Talla-Peso") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))


Anos<-years
SSB<-Sum.Pop.Mod(Pop.Mod,"SSB")
Biomasa=SSB$SSB[,,1]
Reclutamento<-RBH(Biomasa,15000,50)
plot6<-data.frame(Biomasa,Reclutamento)
p6<-ggplot(plot6, aes(x=Biomasa, y=Reclutamento)) +
  geom_line() +
  geom_point() +
  labs(title="Recrutamento",x="Biomasa (t)",y="Recrutamento (t)") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))


plot7<-data.frame(Idades,Mvec)
p7<-ggplot(plot7, aes(x=Idades, y=Mvec)) +
  geom_line() +  coord_cartesian(ylim = c(0, NA)) +
  geom_point() +
  labs(y="Mortalidade natural",title="Mortalidade Natural") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))

p8<-ggplot(plot1, aes(x=Anos, y=Biomasa)) +
  geom_line() +
  labs(x="Anos",y="Biomasa (t)",title="Cortes na Serie temporal da Biomasa") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5)) + 
  geom_vline(xintercept = c(2030,2040,2050),linetype = "dashed",col="red")

prodfun <- function(r,Bt,K,p) return((r*Bt/p)*(1-(Bt/K)^p)) #ec11.2    
densdep <- function(Bt,K,p) return((1/p)*(1-(Bt/K)^p))      
r <- 0.75; K <- 1000.0; Bt <- 1:1000     
sp <- prodfun(r,Bt,K,1.0)  # Schaefer equivalent     
sp0 <- prodfun(r,Bt,K,p=1e-08)  # Fox equivalent
sp0.5<-prodfun(r,Bt,K,0.5)
sp3 <- prodfun(r,Bt,K,3)

sps<-c(sp,sp0,sp0.5,sp3)
type<-c(rep("p=1 Schaefer",length(sp)),rep("p=0 Fox",length(sp)),rep("p=0.5",length(sp)),
        rep("p=3",length(sp)))
Bt<-rep(Bt,4)
plot9<-data.frame(Bt,sps,type)

p9<-ggplot(plot9, aes(x=Bt, y=sps, group=type, color=type)) +
  geom_line() +
  labs(x="Biomasa (t)",y="Produción Excedente (t)",title="Curva de produción") + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))


jpeg("ps.jpeg", width = 2000, height = 1500, res = 300)
p9
dev.off()
