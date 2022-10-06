source("/data/public/apapada1/Scripts/Figures/utilities.R")
library(plyr)


figure2pre<- function(type){
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    return(colnames(slope)[6:ncol(slope)])
}



figure2<-function(type,nx,ty,to,count=0)
{
    if(count==0){count<-type}
    areao<-polII_speed(type,nx,ty,to)
    areao$speedcat <- ifelse(areao$SL < 0, "slower", "faster")
    if(is.numeric(ty)){ty1<-ty}else{
        ty1<-which(grepl(paste0(ty,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[ty1])
    if(is.numeric(to)){to1<-to}else{
        to1<-which(grepl(paste0(to,"."),colnames(areao),fixed=F))}
    print(colnames(areao)[to1])
    if(length(ty1)>1){
        keep<-apply(areao[,ty1],1,min)>apply(areao[,to1],1,max)|apply(areao[,to1],1,min)>apply(areao[,ty1],1,max)
    }else{keep<-areao[,ty1]>apply(areao[,to1],1,max)|apply(areao[,to1],1,min)>areao[,ty1]}
    areao1<-areao[keep,]
    return(list(areao,areao1))
    
    
}
medn<-list()
type="celegans_mut"
figure2pre(type)
nx=c(1:6)
ty="wt"
to="ama1"
medn[[1]] <- figure2(type,nx,ty,to)
names(medn)[[1]]<-"14_ama"

figure2pre(type)
nx=c(10:15)
ty="wtday1"
to="ama1day1"
medn[[2]] <- figure2(type,nx,ty,to,"celegans_mut1")
names(medn)[[2]]<-"1_ama"

figure2pre(type)
nx=c(1:3,10:12)
ty=9:11
to=6:8
medn[[3]] <- figure2(type,nx,ty,to,"celegans_mut2")
names(medn)[[3]]<-"14_wt"
medn<-medn[c(3,1,2)]


slopes<-medn
medn<-sapply(slopes,"[",1)
medn<-lapply(medn,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn,function(x)x$conf.int[[2]]))
library(RColorBrewer)
colors <- brewer.pal(9, "Paired")
colors2 <- brewer.pal(3, "Set3")
ypos <- c(0.4,0.4,0.9)
#pdf(paste0("/data/public/apapada1/Scripts/Figures/Figure2a_unc.p"),width=1000, height=1000)
plot(avg,cex.axis=1.3, ypos, col= c(colors[c(2)], "#F39200",  "#F39200"), xlab="", ylab="", frame.plot=FALSE,  pch=1,  ylim=c(0,1), yaxt = 'n', lwd = 1.2, cex = 2.1, xlim=c(-1.8,1.8))
abline(v=0, lty=2,lwd = 2)
#axis(side = 2, lwd = 1.5, cex = 1.5, labels=FALSE)

text(x = avg-0.2, y = ypos-0.1, labels = c("14 d vs 1 d", "14 d ama-1 vs 14 d", "1 d ama-1 vs 1 d"), pos = 1, xpd = TRUE, cex = 1.2)
mtext("Transcriptional elongation speed(Log2FC)", side = 1, line = 2.2, cex = 1.4)
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)

medn<-sapply(slopes,"[",2)
medn<-lapply(medn,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn,function(x)x$conf.int[[2]]))

points(avg,cex.axis=1.3, ypos, col= c(colors[c(2)], "#F39200",  "#F39200"), xlab="", ylab="", frame.plot=FALSE,  pch=19,  ylim=c(0,1), yaxt = 'n', lwd = 1.2, cex = 2.1, xlim=c(-1.8,1.8))
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)

type="dmelanogaster_mut"
figure2pre(type)
nx=c(1:6)
ty="wt10"
to="wt50"
medn[[1]] <- figure2(type,nx,ty,to)
names(medn)[[1]]<-"50vs10"

figure2pre(type)
nx=c(7:12)
ty="mut10"
to="mut50"
medn[[2]] <- figure2(type,nx,ty,to,"dmelanogaster_mut1")
names(medn)[[2]]<-"50vs10mut"

figure2pre(type)
nx=c(4:6,10:12)
ty="wt50"
to="mut50"
medn[[3]] <- figure2(type,nx,ty,to,"dmelanogaster_mut2")
names(medn)[[3]]<-"50mut"

figure2pre(type)
nx=c(1:3,7:9)
ty="wt10"
to="mut10"
medn[[4]] <- figure2(type,nx,ty,to,"dmelanogaster_mut3")
names(medn)[[4]]<-"10mut"

slopes<-medn
medn<-sapply(slopes,"[",1)
medn<-lapply(medn,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn,function(x)x$conf.int[[2]]))
library(RColorBrewer)
colors <- brewer.pal(9, "Paired")
colors2 <- brewer.pal(3, "Set3")
ypos <- c(0.2,0.6,1,1.35)
plot(avg, ypos, ylim=c(0,1.4),col= c(colors[c(2)],colors[c(2)] ,"#F39200",  "#F39200"), xlab="", ylab="", frame.plot=FALSE,  pch=1, yaxt = 'n', lwd = 1.2, cex = 2.1, xlim=c(-1,1),cex.axis=1.3)
abline(v=0, lty=2,lwd = 2)
#axis(side = 2, lwd = 1.5, cex = 1.5, labels=FALSE)

text(x = avg-0.1, y = ypos, labels = c("50 d vs 10 d\n", "50 d RpII215 vs 10 d RpII215\n", "10 d RpII215 vs 10 d\n", "50 d RpII215 vs 50 d\n"), pos = 1, xpd = TRUE, cex = 1.2)
mtext("Transcriptional elongation speed Log2FC", side = 1, line = 2.2, cex = 1.4)
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
medn<-sapply(slopes,"[",2)
medn<-lapply(medn,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn,function(x)x$conf.int[[2]]))
points(avg, ypos, ylim=c(0,1.4),col= c(colors[c(2)],colors[c(2)] ,"#F39200",  "#F39200"), xlab="", ylab="", frame.plot=FALSE,  pch=19, yaxt = 'n', lwd = 1.2, cex = 2.1, xlim=c(-1,1),cex.axis=1.3)

