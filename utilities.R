speed_filt <- function(type,samplex,sampley, mode="strict")
{
    
    #load(file=paste0("/data/public/apapada1/Revisions/Introns/", type,"_areao"))
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    areao<-areao[,c("ensembl_gene_id",samplex,sampley)]
    areao$speed <- rowMeans(areao[,sampley])/rowMeans(areao[,samplex])
    areao$speedcat <- ifelse(rowMeans(areao[,sampley])/rowMeans(areao[,samplex]) < 1, "slower", "faster")
    keep<-apply(areao[,samplex],1,min)>apply(areao[,sampley],1,max)|apply(areao[,sampley],1,min)>apply(areao[,samplex],1,max)
    areao1<-areao[keep,]
    return(list(areao,areao1))
}

polII_speed <-  function(type, nx=1:6, ty, to, dataset, len=1000)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    # area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    if(length(ty1)>1){
    ynf <- rowMeans(areao[,ty1])}else{ynf=areao[,ty1]}
    if(length(to1)>1){
    onf <- rowMeans(areao[,to1])}else{onf=areao[,to1]}
    areao$young<-ynf
    areao$old<-onf
    areao$SL <- log2(ynf/onf)
    return(areao)
}

polII_scatter <-  function(type, nx=1:6, ty, to,l1,l2, main=0, len=1000)
{
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
  slope <- na.omit(slope)
  
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
  intercept_end <- na.omit(intercept_end)
  
  area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
  samplelength<-ncol(slope)-5
  filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
  # area <-  na.omit(area[filter_reads,])
  ny=nx+samplelength+5
  
  filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
  mean_cov <- apply(area[,nx+5], 2, mean)
  areao <- na.omit(area[filter_reads,])
  areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
  areao<-areao[,c(1:5,nx+5)]
  
  
  areao$length <- areao$end - areao$start
  areao <- areao[areao$length >= len,]
  if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[ty1])
  if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[to1])
  if(length(ty1)>1){
    ynf <- rowMeans(areao[,ty1])}else{ynf=areao[,ty1]}
  if(length(to1)>1){
    onf <- rowMeans(areao[,to1])}else{onf=areao[,to1]}
  areao$young<-log10(-ynf)
  areao$old<-log10(-onf)
  if(main==0){main=type}
  library(ggplot2)
  #pdf(paste0("/data/public/apapada1/Figures/",main,"scatter.pdf"),width=8.5, height=7.5)
  sp <- ggplot(data=areao, aes(x=old, y=young)) + geom_point()+
       geom_abline(intercept = 0, slope = 1,colour='red')+ labs(x=l2,y=l1)
  return(sp)
  #dev.off()
}

polII_speed_corrected <-  function(areao,type)
{
    counts<-read.table(paste0("/data/public/apapada1/Revisions/KallistoNew/featurecounts/",type,"_deseq.tsv"))
    counts<-na.omit(counts)
    colnames(counts)[2]<-"FC"
    print(colnames(counts)[2])
    areao<-merge(areao,counts,by.y=0,by.x=5,all=F)
    areao$SL1<-areao$SL
    areao$SL<-areao$SL+areao$FC
    return(areao)
}


counts_deseq <-  function(type)
{
    counts<-read.table(paste0("/data/public/apapada1/Revisions/KallistoNew/featurecounts/",type,"_deseq.tsv"))
    counts<-na.omit(counts)
    colnames(counts)[2]<-"SL"
    return(counts)
}


polII_speed_pca <-  function(type, nx=1:6, ty, to,main=0,include=T)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    #area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    library(RColorBrewer)
    colors <- brewer.pal(10, "Paired")
    col=c("#4D4D4D",colors[2])
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= 1000,]
    if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    areao.pca<-prcomp(as.matrix(t(areao[,c(ty1,to1)])), center = T,scale. = T)
    if(include==T){
      levs<-as.factor(c(rep("Young",length(ty1)),rep("Old",length(to1))))
      levs<-as.factor(c(rep("Proliferating",length(ty1)),rep("Senescent",length(to1))))
      levs<-relevel(levs,as.character(levs[1]))
      pal=col
    }else{
      #levs<-as.factor(c(rep("WT",length(ty1)),rep("Mutant",length(to1))))
      levs<-as.factor(c(rep("WT",length(ty1)),rep("Progeria",length(to1))))      
      levs<-relevel(levs,as.character(levs[1]))
      pal=c("blue","orange")}
    library(factoextra)
    if(main==0){main<-type}
    pdf(paste0("/data/public/apapada1/Figures/",main,"slopepca.pdf"),width=8.5, height=7.5)
    p<-fviz_pca_ind(areao.pca, geom=c("point"),
                    habillage=levs, addEllipses = F,
                    palette = pal, pointshape=19, labelsize=6,pointsize=10,
                    repel = TRUE ,    # Avoid text overlapping
                    title="",mean.point=F
    ) +
        theme(text = element_text(size = 12),
              plot.title=element_text(size=30), 
              legend.title=element_text(size=40),
              legend.text=element_text(size=30),
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 24))
    print(p)
    dev.off()
   
}


polII_speed_heatmap <-  function(type, nx=1:6, ty, to,main=0,include=T)
{
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
  slope <- na.omit(slope)
  
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
  intercept_end <- na.omit(intercept_end)
  
  area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
  samplelength<-ncol(slope)-5
  filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
  #area <-  na.omit(area[filter_reads,])
  ny=nx+samplelength+5
  
  filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
  mean_cov <- apply(area[,nx+5], 2, mean)
  areao <- na.omit(area[filter_reads,])
  areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
  areao<-areao[,c(1:5,nx+5)]
  
library(pheatmap)
  areao$length <- areao$end - areao$start
  areao <- areao[areao$length >= 1000,]
  if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[ty1])
  if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[to1])
  if(main==0){main<-type}
  areao.heatmap<-as.matrix(areao[,c(ty1,to1)])
  pdf(paste0("/data/public/apapada1/Figures/",main,"slopeheatmap.pdf"),width=8.5, height=7.5)
  areao.cor <- cor(areao.heatmap)
  samples<-gsub(".y","",colnames(areao.cor))
  print(pheatmap(areao.cor,main=main))
  dev.off()
}
slope_exp <-  function(type, nx=1:6, ty, to, dataset, len=1000)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
        ty<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty])
    if(is.numeric(to)){to1<-to}else{
        to<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to])
    if(length(ty)==1){ynf<-areao[,ty]}else{
    ynf <- rowMeans(areao[,ty])}
    if(length(to)==1){onf<-areao[,t0]}else{
    onf <- rowMeans(areao[,to])}
    areao$young<-ynf
    areao$old<-onf
    areao$SL <- log2(ynf/onf)
    return(areao)
}

intron_exp <-  function(type, nx=1:6, ty, to, dataset, len=1000)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/avg.cov.RData"))
    slope <- na.omit(gen_data)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x > 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
        ty<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty])
    if(is.numeric(to)){to1<-to}else{
        to<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to])
    ynf <- rowMeans(areao[,ty])
    onf <- rowMeans(areao[,to])
    areao$young<-ynf
    areao$old<-onf
    areao$SL <- log2(onf/ynf)
    return(areao)
}

intron_slope <-  function(type, nx=1:6, ty, to, dataset, len=1000)
{
    
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/avg.cov.RData"))
    cov <- na.omit(gen_data)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
        ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
        to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    ynf <- rowMeans(areao[,ty1])
    onf <- rowMeans(areao[,to1])
    areao$young<-ynf
    areao$old<-onf
    areao$SL <- log2(ynf/onf)
    if(is.numeric(ty)){ty1<-ty}else{
        ty1<-which(grepl(paste0(ty,"."),colnames(cov),fixed=F))}
    print(colnames(cov)[ty1])
    if(is.numeric(to)){to1<-to}else{
        to1<-which(grepl(paste0(to,"."),colnames(cov),fixed=F))}
    print(colnames(cov)[to1])
    ynf <- rowMeans(cov[,ty1])
    onf <- rowMeans(cov[,to1])
    cov$FC <- log2(onf/ynf)
    areao<-merge(cov,areao,by=c(1:5))
    summary(lm(areao$SL~areao$length))
    summary(lm(log(-areao$YoungR1.y)~log(areao$YoungR1)+log(areao$length)))
  
    

    return(areao)
}

polII_speed_heatmap1 <-  function(type, nx=1:6, ty, to,main=0,include=T)
{
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
  slope <- na.omit(slope)
  
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
  intercept_end <- na.omit(intercept_end)
  
  area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
  samplelength<-ncol(slope)-5
  filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
  #area <-  na.omit(area[filter_reads,])
  ny=nx+samplelength+5
  
  filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
  mean_cov <- apply(area[,nx+5], 2, mean)
  areao <- na.omit(area[filter_reads,])
  areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
  areao<-areao[,c(1:5,nx+5)]
  
  library(pheatmap)
  areao$length <- areao$end - areao$start
  areao <- areao[areao$length >= 1000,]
  if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[ty1])
  if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
  if(main==0){main<-type}
  print(colnames(areao)[to1])
  areao.heatmap<-as.matrix(areao[,c(ty1,to1)])
  pdf(paste0("/data/public/apapada1/Figures/",main,"slopeheatmap1.pdf"),width=8.5, height=7.5)
  areao.cor <- cor(areao.heatmap,method="spearman")
  samples<-gsub(".y","",colnames(areao.cor))
  print(pheatmap(areao.cor,main=main))
  dev.off()
}

polII_speed_heatmap2 <-  function(type, nx=1:6, ty, to,main=0,include=T)
{
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
  slope <- na.omit(slope)
  
  load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
  intercept_end <- na.omit(intercept_end)
  
  area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
  samplelength<-ncol(slope)-5
  filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
  #area <-  na.omit(area[filter_reads,])
  ny=nx+samplelength+5
  
  filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
  mean_cov <- apply(area[,nx+5], 2, mean)
  areao <- na.omit(area[filter_reads,])
  areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
  areao<-areao[,c(1:5,nx+5)]

  
  library(pheatmap)
  areao$length <- areao$end - areao$start
  areao <- areao[areao$length >= 1000,]
  if(is.numeric(ty)){ty1<-ty}else{
    ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
  print(colnames(areao)[ty1])
  if(is.numeric(to)){to1<-to}else{
    to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
  if(main==0){main<-type}
  print(colnames(areao)[to1])
  areao.heatmap<-as.matrix(areao[,c(ty1,to1)])
  colnames(areao.heatmap)<-gsub("Signal.Unique.str2.out.bg","",colnames(areao.heatmap))
   areao.cor <- log10(-areao.heatmap)
  colnames(areao.cor)<-gsub(".y","",colnames(areao.cor))
  return(pheatmap(areao.cor,main="",show_rownames = F,scale="row",treeheight_row=0,fontsize=27))

}

extended2 <-  function(type, nx=1:6, ty, to,l1=ty,l2=to,titl=0, len=1000)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    # area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
        ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
        to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    if(length(ty1)>1){
        ynf <- rowMeans(areao[,ty1])}else{ynf=areao[,ty1]}
    if(length(to1)>1){
        onf <- rowMeans(areao[,to1])}else{onf=areao[,to1]}
    areao$young<-ynf
    areao$old<-onf
    if(titl==0){titl=type}
    #pdf(paste0("/data/public/apapada1/Figures/",titl,"density.pdf"),width=8.5, height=7.5)
    return(plot_density(areao$young,areao$old,l1=l1,l2=l2))
    #dev.off()
    
}

plot_density <-  function(x, y, xlab="density", l1 = "Young", l2= "Old", c1=col[1], c2=col[2], descript="", nb=0, axislab=FALSE)
{
    plot(density(log10(-x), na.rm=T), xlab="", ylab="", cex.axis=1.5 , cex.lab=1.5, frame.plot=FALSE, col=c1, lwd=2, main="")
    #Axis(side=1, labels=axislab)
    lines(density(log10(-y), na.rm=T), col=c2, lwd=2)
    mtext(side=2, xlab, line=2.9,cex=1.5)
    mtext(side=1, "Log10(slope)",line=2.5,cex=1.5)
    legend("topleft", c(l1, l2), lwd=1.2, col=c(c1,c2),  bty = "n", cex=1.5)
    
}


extended2 <-  function(type, nx=1:6, ty, to,l1=ty,l2=to,titl=0, len=1000)
{
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    slope <- na.omit(slope)
    
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/intercept_end.RData"))
    intercept_end <- na.omit(intercept_end)
    
    area <- merge(intercept_end, slope, by=c('chr', 'start', 'end', 'strand', "id"))
    samplelength<-ncol(slope)-5
    filter_reads <- apply(area[,nx+5], 1, function(x)((min(x) > 0)))
    # area <-  na.omit(area[filter_reads,])
    ny=nx+samplelength+5
    
    filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny))))
    mean_cov <- apply(area[,nx+5], 2, mean)
    areao <- na.omit(area[filter_reads,])
    areao<-areao[,-which(grepl(".x",colnames(areao),fixed=T))]
    areao<-areao[,c(1:5,nx+5)]
    
    areao$length <- areao$end - areao$start
    areao <- areao[areao$length >= len,]
    if(is.numeric(ty)){ty1<-ty}else{
        ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
        to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    if(length(ty1)>1){
        ynf <- rowMeans(areao[,ty1])}else{ynf=areao[,ty1]}
    if(length(to1)>1){
        onf <- rowMeans(areao[,to1])}else{onf=areao[,to1]}
    areao$young<-ynf
    areao$old<-onf
    if(titl==0){titl=type}
    #pdf(paste0("/data/public/apapada1/Figures/",titl,"density.pdf"),width=8.5, height=7.5)
    return(plot_density(areao$young,areao$old,l1=l1,l2=l2))
    #dev.off()
    
}

