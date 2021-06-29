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

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/1.-biomasa/B01.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/1.-biomasa/B03.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/1.-biomasa/B05.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3)*100
Clase<-c("0.1","0.1","0.3","0.3","0.5","0.5")
Modelo<-rep(c("ASPIC","SPiCT"),3)
ncs<-data.frame(Proporcion,Clase,Modelo)

p1<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="CV",y="Non-converxentes (%)",title="Erro na observación (Biomasa)") + 
  theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.1.-capturas/B1C1.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.1.-capturas/B1C3.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.1.-capturas/B1C5.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3)*100
Clase<-c("0.1","0.1","0.3","0.3","0.5","0.5")
Modelo<-rep(c("ASPIC","SPiCT"),3)
ncs<-data.frame(Proporcion,Clase,Modelo)

p2<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="CV",y="",title="Erro na observación (Capturas)") +
  theme_minimal() + lims(y=c(NA,15)) +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B1C1.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B3C3.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/2.2.-capturas_B/B5C5.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3)*100
Clase<-c("0.1","0.1","0.3","0.3","0.5","0.5")
Modelo<-rep(c("ASPIC","SPiCT"),3)
ncs<-data.frame(Proporcion,Clase,Modelo)

p3<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="CV's",y="",title="Erro na observación (B + C)") + 
  theme_minimal() + lims(y=c(0,NA)) +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.1.-curta/2001.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.1.-curta/3001.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.1.-curta/4001.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.1.-curta/5001.RData")
nc4<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3,nc4)*100
Clase<-c("20-","20-","30-","30-","40-","40-","50-","50-")
Modelo<-rep(c("ASPIC","SPiCT"),4)
ncs<-data.frame(Proporcion,Clase,Modelo)

p4<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="Período",y="Non-converxentes (%)",title="Series curtas") + 
  theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.2.-curta_B/5001.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.2.-curta_B/5003.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/3.2.-curta_B/5005.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3)*100
Clase<-c("0.1","0.1","0.3","0.3","0.5","0.5")
Modelo<-rep(c("ASPIC","SPiCT"),3)
ncs<-data.frame(Proporcion,Clase,Modelo)

p5<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="CV's",y="",title="Serie curta + Erro observación (Biomasa)") + 
  theme_minimal() +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


load("/Users/anxopaz/Desktop/MTE/TFM/resultados/4.-proceso/P03.RData")
nc1<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/4.-proceso/P05.RData")
nc2<-c(aspic_noconv_prop,spict_noconv_prop)

load("/Users/anxopaz/Desktop/MTE/TFM/resultados/4.-proceso/P07.RData")
nc3<-c(aspic_noconv_prop,spict_noconv_prop)

Proporcion<-c(nc1,nc2,nc3)*100
Clase<-c("0.3","0.3","0.5","0.5","0.7","0.7")
Modelo<-rep(c("ASPIC","SPiCT"),3)
ncs<-data.frame(Proporcion,Clase,Modelo)

p6<-ggplot(ncs, aes(x=Clase, y=Proporcion, group=Modelo, color=Modelo)) + 
  geom_point() +
  geom_line(linetype="dashed") +
  labs(x="CV",y="",title="Erro de Proceso") + 
  theme_minimal() + lims(y=c(0,NA)) +
  theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),
        plot.title=element_text(size=9,face="bold",hjust=0.5))


jpeg("nc.jpeg", width = 2000, height = 1500, res = 300)
grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,ncol=3, nrow=2)
dev.off()

