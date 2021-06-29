library(ggplot2)
library(gridExtra)
library(grid)

grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- switch (position,"bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),legend,ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)),
      "right" = arrangeGrob(do.call(arrangeGrob, gl), legend, ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    invisible(combined)
}

################################################################################

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B1C1.RData")

AA<-matrix(rep(0,11*6),nrow=11,ncol=6)

names<-c("MSY","Bmsy","K","r","n","Bly.Bmsy")
colnames(AA)<-names

value<-c("real","aspic_mean","aspic_lq","aspic_uq","spict_mean","spict_lq","spict_uq",
         "aspic_RMSE","aspic_MAPE","spict_RMSE","spict_MAPE")
rownames(AA)<-value

colnames(spict_estimations)[9]<-("Bly.Bmsy")
colnames(spict_estimations_MAPE)[9]<-("Bly.Bmsy")
colnames(spict_estimations_RMSE)[9]<-("Bly.Bmsy")
colnames(SPM_values)[9]<-("Bly.Bmsy")

for (i in names) {
  AA[1,i]<-SPM_values[,i]
  a_iestimation<-aspic_estimations[,i]
  AA[2,i]<-mean(a_iestimation)
  qq<-quantile(a_iestimation,probs=c(0.05,0.95))
  AA[3,i]<-qq[[1]]
  AA[4,i]<-qq[[2]]
  s_iestimation<-spict_estimations[,i]
  AA[5,i]<-mean(s_iestimation)
  qq<-quantile(s_iestimation,probs=c(0.05,0.95))
  AA[6,i]<-qq[[1]]
  AA[7,i]<-qq[[2]]
  AA[8,i]<-aspic_estimations_RMSE[,i]
  AA[9,i]<-aspic_estimations_MAPE[,i]
  AA[10,i]<-spict_estimations_RMSE[,i]
  AA[11,i]<-spict_estimations_MAPE[,i]
}

AA<-cbind(AA,value)

AA2<-matrix(rep(0,2*12),ncol=2)
names2<-c("B","B/Bmsy")
colnames(AA2)<-names2
values2<-c("aspic_RMSE_mean","aspic_RMSE_lq","aspic_RMSE_uq",
           "aspic_MAPE_mean","aspic_MAPE_lq","aspic_MAPE_uq",
           "spict_RMSE_mean","spict_RMSE_lq","spict_RMSE_uq",
           "spict_MAPE_mean","spict_MAPE_lq","spict_MAPE_uq")
rownames(AA2)<-values2

for (i in names2) {
  AA2[1,i]<-aspic_RMSE["media",i]
  AA2[2,i]<-aspic_RMSE["lq",i]
  AA2[3,i]<-aspic_RMSE["uq",i]
  AA2[4,i]<-aspic_MAPE["media",i]
  AA2[5,i]<-aspic_MAPE["lq",i]
  AA2[6,i]<-aspic_MAPE["uq",i]
  AA2[7,i]<-spict_RMSE["media",i]
  AA2[8,i]<-spict_RMSE["lq",i]
  AA2[9,i]<-spict_RMSE["uq",i]
  AA2[10,i]<-spict_MAPE["media",i]
  AA2[11,i]<-spict_MAPE["lq",i]
  AA2[12,i]<-spict_MAPE["uq",i]
}


################################################################################

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B3C3.RData")

BB<-matrix(rep(0,11*6),nrow=11,ncol=6)
colnames(BB)<-names
rownames(BB)<-value

colnames(spict_estimations)[9]<-("Bly.Bmsy")
colnames(spict_estimations_MAPE)[9]<-("Bly.Bmsy")
colnames(spict_estimations_RMSE)[9]<-("Bly.Bmsy")
colnames(SPM_values)[9]<-("Bly.Bmsy")

for (i in names) {
  BB[1,i]<-SPM_values[,i]
  a_iestimation<-aspic_estimations[,i]
  BB[2,i]<-mean(a_iestimation)
  qq<-quantile(a_iestimation,probs=c(0.05,0.95))
  BB[3,i]<-qq[[1]]
  BB[4,i]<-qq[[2]]
  s_iestimation<-spict_estimations[,i]
  BB[5,i]<-mean(s_iestimation)
  qq<-quantile(s_iestimation,probs=c(0.05,0.95))
  BB[6,i]<-qq[[1]]
  BB[7,i]<-qq[[2]]
  BB[8,i]<-aspic_estimations_RMSE[,i]
  BB[9,i]<-aspic_estimations_MAPE[,i]
  BB[10,i]<-spict_estimations_RMSE[,i]
  BB[11,i]<-spict_estimations_MAPE[,i]
}

BB<-cbind(BB,value)


BB2<-matrix(rep(0,2*12),ncol=2)
colnames(BB2)<-names2
rownames(BB2)<-values2

for (i in names2) {
  BB2[1,i]<-aspic_RMSE["media",i]
  BB2[2,i]<-aspic_RMSE["lq",i]
  BB2[3,i]<-aspic_RMSE["uq",i]
  BB2[4,i]<-aspic_MAPE["media",i]
  BB2[5,i]<-aspic_MAPE["lq",i]
  BB2[6,i]<-aspic_MAPE["uq",i]
  BB2[7,i]<-spict_RMSE["media",i]
  BB2[8,i]<-spict_RMSE["lq",i]
  BB2[9,i]<-spict_RMSE["uq",i]
  BB2[10,i]<-spict_MAPE["media",i]
  BB2[11,i]<-spict_MAPE["lq",i]
  BB2[12,i]<-spict_MAPE["uq",i]
}


################################################################################

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B5C5.RData")

CC<-matrix(rep(0,11*6),nrow=11,ncol=6)
colnames(CC)<-names
rownames(CC)<-value

colnames(spict_estimations)[9]<-("Bly.Bmsy")
colnames(spict_estimations_MAPE)[9]<-("Bly.Bmsy")
colnames(spict_estimations_RMSE)[9]<-("Bly.Bmsy")
colnames(SPM_values)[9]<-("Bly.Bmsy")

for (i in names) {
  CC[1,i]<-SPM_values[,i]
  a_iestimation<-aspic_estimations[,i]
  CC[2,i]<-mean(a_iestimation)
  qq<-quantile(a_iestimation,probs=c(0.05,0.95))
  CC[3,i]<-qq[[1]]
  CC[4,i]<-qq[[2]]
  s_iestimation<-spict_estimations[,i]
  CC[5,i]<-mean(s_iestimation)
  qq<-quantile(s_iestimation,probs=c(0.05,0.95))
  CC[6,i]<-qq[[1]]
  CC[7,i]<-qq[[2]]
  CC[8,i]<-aspic_estimations_RMSE[,i]
  CC[9,i]<-aspic_estimations_MAPE[,i]
  CC[10,i]<-spict_estimations_RMSE[,i]
  CC[11,i]<-spict_estimations_MAPE[,i]
}

CC<-cbind(CC,value)


CC2<-matrix(rep(0,2*12),ncol=2)
colnames(CC2)<-names2
rownames(CC2)<-values2

for (i in names2) {
  CC2[1,i]<-aspic_RMSE["media",i]
  CC2[2,i]<-aspic_RMSE["lq",i]
  CC2[3,i]<-aspic_RMSE["uq",i]
  CC2[4,i]<-aspic_MAPE["media",i]
  CC2[5,i]<-aspic_MAPE["lq",i]
  CC2[6,i]<-aspic_MAPE["uq",i]
  CC2[7,i]<-spict_RMSE["media",i]
  CC2[8,i]<-spict_RMSE["lq",i]
  CC2[9,i]<-spict_RMSE["uq",i]
  CC2[10,i]<-spict_MAPE["media",i]
  CC2[11,i]<-spict_MAPE["lq",i]
  CC2[12,i]<-spict_MAPE["uq",i]
}


################################################################################
################################################################################
################################################################################

plotting<-rbind(AA,BB,CC)
plotting2<-rbind(AA2,BB2,CC2)
colnames(plotting)[6]<-"BT/Bmsy"
names[6]<-"BT/Bmsy"

clase<-c(rep("0.1",length(value)),
  rep("0.3",length(value)),
  rep("0.5",length(value)))
plotting<-cbind(plotting,clase)

clase2<-c(rep("0.1",length(values2)),
         rep("0.3",length(values2)),
         rep("0.5",length(values2)))
value2<-rep(values2,3)
plotting2<-cbind(plotting2,value2,clase2)


plo<-list()
for (i in names){
  inde<-which(plotting[,7]=="real"|plotting[,7]=="aspic_mean"|plotting[,7]=="spict_mean")
  medias<-as.numeric(plotting[inde,i])
  inde2<-which(plotting[,7]=="real"|plotting[,7]=="aspic_lq"|plotting[,7]=="spict_lq")
  lq<-as.numeric(plotting[inde2,i])
  cvs<-factor(c(rep(0.1,3),rep(0.3,3),rep(0.5,3)))
  valores<-rep(c("real","ASPIC","SPiCT"),3)
  ceros<-which(valores=="real")
  lq[ceros]=0
  inde3<-which(plotting[,7]=="real"|plotting[,7]=="aspic_uq"|plotting[,7]=="spict_uq")
  uq<-as.numeric(plotting[inde3,i])
  cvs<-factor(c(rep(0.1,3),rep(0.3,3),rep(0.5,3)))
  valores<-rep(c("real","ASPIC","SPiCT"),3)
  ceros<-which(valores=="real")
  uq[ceros]=0
  plo[[i]]<-data.frame(cvs=cvs,valores=valores,medias=medias,lq=lq,uq=uq)
}

plo2<-list()
for (i in names){
  inde<-which(plotting[,7]=="aspic_RMSE"|plotting[,7]=="spict_RMSE")
  RMSE<-as.numeric(plotting[inde,i])
  inde2<-which(plotting[,7]=="aspic_MAPE"|plotting[,7]=="spict_MAPE")
  MAPE<-as.numeric(plotting[inde2,i])
  cvs<-factor(c(rep(0.1,2),rep(0.3,2),rep(0.5,2)))
  valores<-rep(c("ASPIC","SPiCT"),3)
  plo2[[i]]<-data.frame(cvs=cvs,valores=valores,RMSE=RMSE,MAPE=MAPE)
}

plo3<-list()
for (i in names2){
  inde<-which(plotting2[,3]=="aspic_RMSE_mean"|plotting2[,3]=="spict_RMSE_mean")
  RMSE_mean<-as.numeric(plotting2[inde,i])
  inde2<-which(plotting2[,3]=="aspic_RMSE_lq"|plotting2[,3]=="spict_RMSE_lq")
  RMSE_lq<-as.numeric(plotting2[inde2,i])
  inde22<-which(plotting2[,3]=="aspic_RMSE_uq"|plotting2[,3]=="spict_RMSE_uq")
  RMSE_uq<-as.numeric(plotting2[inde22,i])
  inde3<-which(plotting2[,3]=="aspic_MAPE_mean"|plotting2[,3]=="spict_MAPE_mean")
  MAPE_mean<-as.numeric(plotting2[inde3,i])
  inde4<-which(plotting2[,3]=="aspic_MAPE_lq"|plotting2[,3]=="spict_MAPE_lq")
  MAPE_lq<-as.numeric(plotting2[inde4,i])
  inde44<-which(plotting2[,3]=="aspic_MAPE_uq"|plotting2[,3]=="spict_MAPE_uq")
  MAPE_uq<-as.numeric(plotting2[inde44,i])
  cvs<-factor(c(rep(0.1,2),rep(0.3,2),rep(0.5,2)))
  valores<-rep(c("ASPIC","SPiCT"),3)
  plo3[[i]]<-data.frame(cvs=cvs,valores=valores,RMSE_mean=RMSE_mean,
                        RMSE_lq=RMSE_lq,RMSE_uq=RMSE_uq,
                        MAPE_mean=MAPE_mean,MAPE_lq=MAPE_lq,MAPE_uq=MAPE_uq)
}

plot_list<-list()
for(i in 1:6) {
  p<-ggplot(plo[[i]], aes(x=cvs, y=medias, group=valores, color=valores)) + 
  geom_point(position=position_dodge(0.20)) + 
  geom_errorbar(data=plo[[i]][-c(1,4,7),],aes(x=cvs,group=valores, color=valores,
                                              ymin=lq, ymax=uq), 
                inherit.aes = FALSE,
                width=.2, position=position_dodge(0.25)) +
  labs(x=element_blank(),y=names(plo)[[i]]) + theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"))
        # legend.background = element_rect(fill="grey100",
        #                                  size=0.1, linetype="solid"))
  plot_list<-c(plot_list,list(p))
}

for(i in c(5,6)) { 
  plot_list[[i]]<-plot_list[[i]] + labs(x="CV's")
}

jpeg("capbio1.jpeg", width = 2000, height = 1500, res = 300)
grid_arrange_shared_legend(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                           plot_list[[4]],plot_list[[5]],plot_list[[6]],
                           ncol=2, nrow=3)
dev.off()


plot_list2<-list()
for(i in 1:6) {
  p<-ggplot(plo2[[i]], aes(x=cvs, y=RMSE, group=valores, color=valores)) + 
    geom_point(position=position_dodge(0)) +
    geom_line(linetype="dashed") +
    labs(x=element_blank(),y=names(plo2)[[i]]) + theme_minimal() +
    theme(legend.title = element_blank(),axis.line = element_line(colour = "black"))
  plot_list2<-c(plot_list2,list(p))
}

plot_list3<-list()
for(i in 1:6) {
  p<-ggplot(plo2[[i]], aes(x=cvs, y=MAPE, group=valores, color=valores)) + 
    geom_point(position=position_dodge(0)) +
    geom_line(linetype="dashed") +
    labs(x=element_blank(),y=names(plo2)[[i]]) + theme_minimal() +
    theme(legend.title = element_blank(),axis.line = element_line(colour = "black"))
  plot_list3<-c(plot_list3,list(p))
}

for(i in c(5,6)) { 
  plot_list2[[i]]<-plot_list2[[i]] + labs(x="CV's")
  plot_list3[[i]]<-plot_list3[[i]] + labs(x="CV's")
}

for(i in c(1,2)) { 
  plot_list2[[i]]<-plot_list2[[i]] + labs(title="RMSE")
  plot_list3[[i]]<-plot_list3[[i]] + labs(title="MAPE")
}

jpeg("capbio2.jpeg", width = 2000, height = 1500, res = 300)
grid_arrange_shared_legend(plot_list2[[1]],plot_list3[[1]],plot_list2[[2]],
                           plot_list3[[2]],plot_list2[[3]],plot_list3[[3]],
                           plot_list2[[4]],plot_list3[[4]],plot_list2[[5]],
                           plot_list3[[5]],plot_list2[[6]],plot_list3[[6]],
                           ncol=4, nrow=3)
dev.off()

plot_list4<-list()
for(i in 1:2) {
  p<-ggplot(plo3[[i]], aes(x=cvs, y=RMSE_mean, group=valores, color=valores)) + 
    geom_point(position=position_dodge(0.25)) +
    geom_errorbar(aes(x=cvs,group=valores, color=valores,
                                                ymin=RMSE_lq, 
                                                ymax=RMSE_uq),
                  width=.2, position=position_dodge(0.25)) +
    geom_line(linetype="dashed") +
    labs(x=element_blank(),y=names(plo3)[[i]]) + theme_minimal() +
    theme(legend.title = element_blank(),axis.line = element_line(colour = "black"))
  plot_list4<-c(plot_list4,list(p))
}

plot_list5<-list()
for(i in 1:2) {
  p<-ggplot(plo3[[i]], aes(x=cvs, y=MAPE_mean, group=valores, color=valores)) + 
    geom_point(position=position_dodge(0.25)) +
    geom_errorbar(aes(x=cvs,group=valores, color=valores,
                      ymin=MAPE_lq, 
                      ymax=MAPE_uq),
                  width=.2, position=position_dodge(0.25)) +
  geom_line(linetype="dashed") +
    labs(x=element_blank(),y=names(plo3)[[i]]) + theme_minimal() +
    theme(legend.title = element_blank(),axis.line = element_line(colour = "black"))
  plot_list5<-c(plot_list5,list(p))
}

for(i in 2) { 
  plot_list4[[i]]<-plot_list4[[i]] + labs(x="CV's")
  plot_list5[[i]]<-plot_list5[[i]] + labs(x="CV's")
}

for(i in 1) { 
  plot_list4[[i]]<-plot_list4[[i]] + labs(title="RMSE")
  plot_list5[[i]]<-plot_list5[[i]] + labs(title="MAPE")
}

jpeg("capbio3.jpeg", width = 2000, height = 1500, res = 300)
grid_arrange_shared_legend(plot_list4[[1]],plot_list5[[1]],plot_list4[[2]],
                           plot_list5[[2]],
                           ncol=2, nrow=2)
dev.off()
