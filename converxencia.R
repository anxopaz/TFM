load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.2.-curta_B/5005.RData")

plots<-function(dat){
  desc=""; fname="convPlot"; step=1; upPerc=0.95; lowPerc=0.05
  dat <- as.numeric(dat)
  datLen <- length(dat)
  datGrid <- seq(step, datLen, by=step)
  cumMean <- unlist(lapply(datGrid, function (x) mean(dat[1:x])))
  cumMedian <- unlist(lapply(datGrid, function (x) median(dat[1:x])))
  cumSD <- unlist(lapply(datGrid, function (x) sd(dat[1:x])))
  cumCV <- cumSD/cumMean
  cumUpper <-  unlist(lapply(datGrid, function (x) quantile(dat[1:x], probs=upPerc)))
  cumLow <-  unlist(lapply(datGrid, function (x) quantile(dat[1:x], probs=lowPerc)))
  cumDF <- data.frame(cumMean,cumSD, cumCV, cumUpper, cumMedian, cumLow)
  names(cumDF) <- c(" media", " sd", " CV", paste("perc ", as.character(upPerc), sep=""),
                    " mediana", paste("perc ", as.character(lowPerc), sep=""))
  
  ## Plots
  par(mfcol=c(1,3),oma=c(0,0,2,0) )
  
  #MEAN
  y=cumDF[,1]
  conv <- y[length(y)]
  plot(datGrid, y, xlab="iters", ylab="media", type='l', lwd=2,
       ylim=c(conv*0.95, conv*1.05))
  abline(h=conv*0.975, lty=1, lwd=1, col=2)
  abline(h=conv*1.025, lty=1, lwd=1, col=2)
  text(x=datLen, y=conv*1.03, labels="5%", cex=0.6)
  abline(h=conv*0.995, lty=2, lwd=.5, col=2)
  abline(h=conv*1.005, lty=2, lwd=.5, col=2)
  text(x=datLen, y=conv*1.008, labels="1%", cex=0.6)
  
  # perc 0.95
  y=cumDF[,4]
  conv <- y[length(y)]
  plot(datGrid, y, xlab="iters", ylab="perc 0.95", type='l', lwd=2,
       ylim=c(conv*0.95, conv*1.05))
  abline(h=conv*0.975, lty=1, lwd=1, col=2)
  abline(h=conv*1.025, lty=1, lwd=1, col=2)
  text(x=datLen, y=conv*1.03, labels="5%", cex=0.6)
  abline(h=conv*0.995, lty=2, lwd=.5, col=2)
  abline(h=conv*1.005, lty=2, lwd=.5, col=2)
  text(x=datLen, y=conv*1.008, labels="1%", cex=0.6)
  
  #perc 0.05
  y=cumDF[,6]
  conv <- y[length(y)]
  plot(datGrid, y, xlab="iters", ylab="perc 0.05", type='l', lwd=2,
       ylim=c(conv*0.95, conv*1.05))
  abline(h=conv*0.975, lty=1, lwd=1, col=2)
  abline(h=conv*1.025, lty=1, lwd=1, col=2)
  text(x=datLen, y=conv*1.03, labels="5%", cex=0.6)
  abline(h=conv*0.995, lty=2, lwd=.5, col=2)
  abline(h=conv*1.005, lty=2, lwd=.5, col=2)
  text(x=datLen, y=conv*1.008, labels="1%", cex=0.6)
}



###### ASPIC ######

aMAPE<-matrix(rep(0,nrep),ncol=2,nrow=nrep)
aRMSE<-matrix(rep(0,nrep),ncol=2,nrow=nrep)
for(i in 1:nrep){
  predicted<-aspic_states$B0est[i,]
  real<-B_SPM
  aMAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$BBmsy[i,]
  real<-B_SPM/Bmsy
  aMAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-aspic_states$B0est[i,]
  real<-B_SPM
  aRMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-aspic_states$BBmsy[i,]
  real<-B_SPM/Bmsy
  aRMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
}


nomes<-c("MSY","Bmsy","K","n","r","BT/Bmsy")
a_estimacions<-aspic_estimations[,c(1,3,4,7,6,9)]
colnames(a_estimacions)<-nomes

jpeg("c5b510.jpeg", width = 2000, height = 600, res = 300)
plots(aMAPE[,2])
dev.off()

for (i in nomes){
  cat(i)
  dat=a_estimacions[,i]
  plots(dat)
}

cat("aMAPE_B")
plots(dat=aMAPE[,1])

cat("aMAPE_BR")
plots(dat=aMAPE[,2])

cat("aRMSE_B")
plots(dat=aRMSE[,1])

cat("aRMSE_BR")
plots(dat=aRMSE[,2])



###### SPiCT ######

sMAPE<-matrix(rep(0,nrep),ncol=2,nrow=nrep)
sRMSE<-matrix(rep(0,nrep),ncol=2,nrow=nrep)
for(i in 1:nrep){
  predicted<-spict_states$B[i,]
  real<-B_SPM
  sMAPE[i,1]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$BR[i,]
  real<-B_SPM/Bmsy
  sMAPE[i,2]<-1/length(years)*(sum(abs((predicted-real)/real)))
  predicted<-spict_states$B[i,]
  real<-B_SPM
  sRMSE[i,1]<-sqrt(sum((predicted-real)^2)/length(years)) 
  predicted<-spict_states$BR[i,]
  real<-B_SPM/Bmsy
  sRMSE[i,2]<-sqrt(sum((predicted-real)^2)/length(years)) 
}


nomes<-c("MSY","Bmsy","K","n","r","BT/Bmsy")
s_estimacions<-spict_estimations[,c(1,3,4,7,6,9)]
colnames(s_estimacions)<-nomes

jpeg("cb5msy.jpeg", width = 2000, height = 600, res = 300)
plots(s_estimacions[,"MSY"])
dev.off()

for (i in 1:6){
  cat(i)
  dat=s_estimacions[,i]
  jpeg(filename=paste0("c5b51",i,".jpeg",sep=""), width = 2000, height = 600, res = 300)
  plots(dat)
  dev.off()
}

cat("sMAPE_B")
plots(dat=sMAPE[,1])

cat("sMAPE_BR")
plots(dat=sMAPE[,2])

cat("sRMSE_B")
plots(dat=sRMSE[,1])

cat("sRMSE_BR")
plots(dat=sRMSE[,2])


###### saves ######

save(aMAPE,aRMSE,sMAPE,sRMSE,a_estimacions,s_estimacions,nomes,plots,
     file="~/Desktop/MTE/TFM/resultados/converxencia/data/5005")

jpeg("c5b520.jpeg", width = 2000, height = 600, res = 300)
plots(sMAPE[,2])
dev.off()
