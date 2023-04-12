setwd("/data/public/apapada1/Revisions/h3")
nx=paste0(c("H3_r1","H3_r2","wt_r1","wt_r2"),".x")
ny=paste0(c("H3_r1","H3_r2","wt_r1","wt_r2"),".y")
ty<-c("wt_r1.y","wt_r2.y")
to<-c("H3_r1.y","H3_r2.y")
nx=paste0(c("H4_r1","H4_r2","wt_r1","wt_r2"),".x")
ny=paste0(c("H4_r1","H4_r2","wt_r1","wt_r2"),".y")
ty<-c("wt_r1.y","wt_r2.y")
to<-c("H4_r1.y","H4_r2.y")
nx=paste0(c("H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2"),".x")
ny=paste0(c("H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2"),".y")

load("gen_diff_slope")
colnames(area_diff)<-c('chr', 'start', 'end', 'strand', 'ensembl_gene_id',"H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2")
area<-area_diff

load("intercept")

colnames(area_diff)<-c('chr', 'start', 'end', 'strand', 'ensembl_gene_id',"H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2")
head(area_diff)


min_cov <- area_diff    
area <- merge(min_cov, area, by=c('chr', 'start', 'end', 'strand','ensembl_gene_id'))
filter_reads <- apply(area[,nx], 1, function(x)((min(x) > 0)))
#area <-  na.omit(area[filter_reads,])
filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
areao <- na.omit(area[filter_reads,])
#load(file="/data/public/cdebes/workspace/transcription_aging/premature/slope_model")
#areao[,ny] <- apply(areao[,ny], 2, function(row)slope_model$coefficients[2]*-log(abs(row)) + slope_model$coefficients[1])
areao$length <- areao$end - areao$start
ynf <- rowMeans(areao[,ty])
to<-c("H3_r1.y","H3_r2.y")
onf <- rowMeans(areao[,to])
areao$control<-ynf
areao$H3<-onf
to<-c("H4_r1.y","H4_r2.y")
onf <- rowMeans(areao[,to])
areao$H4<-onf
areao$SL <- onf/ynf
wilcox.test(onf,ynf,paired=T,conf.int=T)
wilcox.test(onf,areao$wt_r1.y,paired=T,conf.int=T)
wilcox.test(onf,areao$wt_r2.y,paired=T,conf.int=T)
wilcox.test(areao$H3_r1.y,areao$wt_r1.y,paired=T,conf.int=T)
wilcox.test(areao$H3_r1.y,areao$wt_r2.y,paired=T,conf.int=T)
wilcox.test(areao$H3_r2.y,areao$wt_r1.y,paired=T,conf.int=T)
wilcox.test(areao$H3_r2.y,areao$wt_r2.y,paired=T,conf.int=T)
wilcox.test(areao$H4_r1.y,areao$wt_r1.y,paired=T,conf.int=T)
wilcox.test(areao$H4_r1.y,areao$wt_r2.y,paired=T,conf.int=T)
wilcox.test(areao$H4_r2.y,areao$wt_r1.y,paired=T,conf.int=T)
wilcox.test(areao$H4_r2.y,areao$wt_r2.y,paired=T,conf.int=T)

library(factoextra)
pcaplot<-function(type,ty,to,ny,no,plotname,title1,include=F,slopes=T){
  
  areao<-subset(areao,length>1000)
  #  if(include==T){
  # areao.pca<-prcomp(as.matrix(t(areao[,c(ty,to,ny,no)])), center = T,scale. = T)
  #  levs<-as.factor(c(rep(ny,length(ty)),rep(no,length(to)),ny,no))
  # 
  #  }else{
  areao.pca<-prcomp(as.matrix(t(areao[,c(ty,to)])), center = T,scale. = T)
  levs<-as.factor(c(rep(ny,length(ty)),rep(no,length(to))))
  levs<-relevel(levs,ny)
  if(include==T){
    pal=c("blue","orange")
  }else{
    pal=c("red","blue")}
  p<-fviz_pca_ind(areao.pca, geom=c("point"),
                  habillage=levs, addEllipses = F,
                  palette = pal, pointshape=19, labelsize=6,pointsize=10,
                  repel = TRUE ,    # Avoid text overlapping
                  title=title1,mean.point=F
  ) +
    theme(text = element_text(size = 12),
          plot.title=element_text(size=30), 
          legend.title=element_text(size=40),
          legend.text=element_text(size=30),
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 24))
  
  return(p)
  
}
ty<-c("H3_r1.y","H3_r2.y")
to<-c("H3_r1.y","H3_r2.y")
ty<-c("wt_r1.y","wt_r2.y")
to<-c("H4_r1.y","H4_r2.y")
tz<-c("H3_r1.y","H3_r2.y")
pcaplot("mmusculus",ty,to,"Control","H3","mm_pca1.png","H3 overexpression",T,T)
pcaplot("mmusculus",ty,to,"Control","H4","mm_pca1.png","H4 overexpression",T,T)
pcaplot("mmusculus",ty,to,"H3","H4","mm_pca1.png","H4 overexpression",T,T)
pcaplot("mmusculus",ty,to,tz,"Control","H3","H4","mm_pca1.png","Histone overexpression",T,T)
pcaplot<-function(type,ty,to,tz,ny,no,nz,plotname,title1,include=F,slopes=T){
  
  areao<-subset(areao,length>1000)
  #  if(include==T){
  # areao.pca<-prcomp(as.matrix(t(areao[,c(ty,to,ny,no)])), center = T,scale. = T)
  #  levs<-as.factor(c(rep(ny,length(ty)),rep(no,length(to)),ny,no))
  # 
  #  }else{
  areao.pca<-prcomp(as.matrix(t(areao[,c(ty,to,tz)])), center = T,scale. = T)
  levs<-as.factor(c(rep(ny,length(ty)),rep(no,length(to)),rep(nz,length(tz))))
  levs<-relevel(levs,ny)
  pal=c("blue","orange","gray")
  p<-fviz_pca_ind(areao.pca, geom=c("point"),
                  habillage=levs, addEllipses = F,
                  palette = pal, pointshape=19, labelsize=6,pointsize=10,
                  repel = TRUE ,    # Avoid text overlapping
                  title=title1,mean.point=F
  ) +
    theme(text = element_text(size = 12),
          plot.title=element_text(size=30), 
          legend.title=element_text(size=40),
          legend.text=element_text(size=30),
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 24))
  
  return(p)
  
}

areaoh3<-areao
areaoh4<-areao

#Old style
library(RColorBrewer)
colors <- brewer.pal(4, "Paired")
#devSVG("data/public/apapada1/rna_index.svg", width=5.0472, height=10.4)
#par(mfrow=c(1,2))
ypos <- c(0.2, 0.4)
h3<-wilcox.test(areaoh3$H3,areaoh3$control,paired=T,conf.int=T)
h4<-wilcox.test(areaoh4$H4,areaoh4$control,paired=T,conf.int=T)
h3<-wilcox.test(log2((1/areao$H3)/(1/areao$control)),conf.int=T)
h4<-wilcox.test(log2((1/areao$H4)/(1/areao$control)),conf.int=T)
minus <- c(h3$conf.int[[1]], h4$conf.int[[1]])
plus <- c(h3$conf.int[[2]], h4$conf.int[[2]])
avg <- c(h3$estimate[[1]], h4$estimate[[1]])

plot(avg, ypos, xlab="", ylab="", yaxt = 'n',ylim=c(0,0.6) ,xlim=c(-0.8,0.2),pch=19, col= c("orange","grey"),frame.plot=FALSE, cex.axis=1.2,cex=2)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
#axis(side = 2, lwd = 1, labels=FALSE)
mtext("Log2FC of elongation rate", side = 1, line = 2.0, cex = 1.2)
abline(v=0, lwd=2, lty=2)
text(y = ypos, x = c((avg[1]-0.25),avg[2]), labels =  c( "H3-GFP overexpression","H4-GFP overexpression") ,pos = 1, xpd = TRUE, cex = 1.1)
dev.off()

boxplot(areao$control,areao$H3,areao$H4,names=c("Control","H3","H4"),outline=F)
boxplot(log10(-areao$control),log10(-areao$H3),log10(-areao$H4),names=c("Control","H3","H4"),outline=F)

#New style

library(RColorBrewer)
colors <- brewer.pal(4, "Paired")

ypos <- c(0.4, 0.6)
areao<-areao1
areao<-areao[-c(36),]
h3<-log2((1/areao$H3)/(1/areao$wt_r1.y))
h4<-log2((1/areao$H4)/(1/areao$wt_r1.y))
h3<-log2((1/areao$H3)/(1/areao$wt_r2.y))
h4<-log2((1/areao$H4)/(1/areao$wt_r2.y))
h3<-log2((1/areao$H3)/(1/areao$control))
h4<-log2((1/areao$H4)/(1/areao$control))
h3<-t.test(h3,mu=0)
h4<-t.test(h4,mu=0)
h3<-wilcox.test(h3,mu=0,conf.int=T)
h4<-wilcox.test(h4,mu=0,conf.int=T)
h3

minus <- c(h3$conf.int[[1]], h4$conf.int[[1]])
plus <- c(h3$conf.int[[2]], h4$conf.int[[2]])
avg <- c(h3$estimate[[1]], h4$estimate[[1]])

pdf("h3h4.pdf")
plot(avg, ypos, xlab="", ylab="", yaxt = 'n',ylim=c(0.2,0.6) ,xlim=c(-1,1),pch=19, col= c("orange"),frame.plot=FALSE, cex.axis=1.4)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
#axis(side = 2, lwd = 1, labels=FALSE)
mtext("Log2 fold change of transcriptional speed", side = 1, line = 2.0, cex = 1.2)
abline(v=0, lwd=2, lty=2)
text(y = ypos, x = avg, labels =  c( "H3 overexpression","H4 overexpression") ,pos = 1, xpd = TRUE, cex = 1.1)
dev.off()

boxplot(areao$control,areao$H3,areao$H4,names=c("Control","H3","H4"),outline=F)
boxplot(log10(-areao$control),log10(-areao$H3),log10(-areao$H4),names=c("Control","H3","H4"),outline=F)


test1<-log2((1/areao$H3)/(1/areao$control))
test2<-log2((1/areao$H3)/(1/areao$wt_r1.y))
test3<-log2((1/areao$H3)/(1/areao$wt_r2.y))
test<-cbind.data.frame(test1,test2,test3)
colnames(test)<-c("Mean","Wt1","Wt2")

h3<-log2((1/areao$H3)/(1/areao$control))
log2((1/areao$H3)/(1/areao$control))==log2((1/areao$H3)/(1/areao$control))
