#Loading functions
source("/data/public/apapada1/Scripts/Figures/utilities.R")
library(plyr)

#Check if names are right
figure1pre<- function(type){
    load(file=paste0("/data/public/apapada1/Revisions/KallistoNew/",type,"/slope.RData"))
    return(colnames(slope)[6:ncol(slope)])
}

#Filter for consistent introns
extended3<-function(type,nx,ty,to,count=0)
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
        keep<-apply(areao[,ty1],1,min)>apply(areao[,to1],1,max)|apply(areao[,to1],1,min)>apply(areao[,ty1],1,max)
        areao1<-areao[keep,]
        return(list(areao,areao1))
  
    
}

#Calculating slopes by sample
slopes<-list()
type="celegans_mut"
nx=c(1:6)
ty="wt"
to="ama1"
slopes[[1]] <- extended3(type,nx,ty,to)
names(slopes)[[1]]<-"14_ama"

nx=c(10:15)
ty="wtday1"
to="ama1day1"
slopes[[2]] <- extended3(type,nx,ty,to,"celegans_mut1")
names(slopes)[[2]]<-"1_ama"

type="hsapiens_nascent"
nx=1:4
ty="Young"
to="Old"
slopes[[3]] <- extended3(type,nx,ty,to)
names(slopes)[[3]]<-"IMR90"

type="hsapiens"
nx=1:6
ty="young"
to="old"
slopes[[4]] <- extended3(type,nx,ty,to)
names(slopes)[[4]]<-"HUVEC"

type="mmusculus_blood"
nx=1:6
ty=6:8
to=9:11
slopes[[5]] <- extended3(type,nx,ty,to)
#slopes[[5]]<-NA
names(slopes)[[5]]<-"Mouse(blood)"

type="hsapiens_blood"

nx=1:12
ty="Young"
to="Old"
slopes[[6]] <- extended3(type,nx,ty,to)
names(slopes)[[6]]<-"Human(blood)"

type="hsapiens_prog"

nx=1:4
ty="wt"
to="PROG"
slopes[[7]] <- extended3(type,nx,ty,to)
names(slopes)[[7]]<-"Human(Progeria)"

type="mmusculus"

nx=1:6
ty="young"
to="old"
slopes[[8]] <- extended3(type,nx,ty,to)
names(slopes)[[8]]<-"Mouse(kidney)"

type="mmusculus_kidney_CR"

nx=1:8
ty="C2"
to="CR"
slopes[[9]] <- extended3(type,nx,ty,to)
names(slopes)[[9]]<-"Mouse(kidney, CR)"

type="mmusculus_hypothalamus"

nx=1:7
ty="wt"
to="mut"
slopes[[10]] <- extended3(type,nx,ty,to)
names(slopes)[[10]]<-"Mouse(hypothalamus)"

type="rnorvegicus"

nx=1:6
ty="Young"
to="Old"
slopes[[11]] <- extended3(type,nx,ty,to)
names(slopes)[[11]]<-"Rat"

type="celegans"

nx=c(1,4,5,8,9,12)
ty="0Signal"
to="21Signal"
slopes[[12]] <- extended3(type,nx,ty,to)
names(slopes)[[12]]<-type

nx=c(1,3,5,7,9,11)
ty="0Signal"
to="14Signal"
slopes[[13]] <- extended3(type,nx,ty,to)
names(slopes)[[13]]<-"celegans_14"

type = "dmelanogaster_merged"

nx=c(1,4,6,7,11,12)
ty="wt10"
to="wt50"
slopes[[14]] <- extended3(type,nx,ty,to)
names(slopes)[[14]]<-type

nx=c(1,2,3,5,11,12)
ty="wt50"
to="mut50"
slopes[[15]] <- extended3(type,nx,ty,to)
names(slopes)[[15]]<-"dmelanogaster_dilp"

type="mmusculus_aldr"

nx=1:6
ty="AL5"
to="DR5"
slopes[[16]] <- extended3(type,nx,ty,to)
names(slopes)[[16]]<-"ALDR5"

type="mmusculus_aldr"

nx=7:12
ty="AL16"
to="DR16"
slopes[[17]] <- extended3(type,nx,ty,to)
names(slopes)[[17]]<-"ALDR16"

type="mmusculus_aldr"

nx=13:18
ty="AL26.5"
to="DR26.5"
slopes[[18]] <- extended3(type,nx,ty,to)
names(slopes)[[18]]<-type

type="celegans_mut"

nx=c(1:3,7:9)
ty="wt"
to="daf2"
slopes[[19]] <- extended3(type,nx,ty,to,"celegans_daf2")
names(slopes)[[19]]<-"daf2"

type = "dmelanogaster_merged"

nx=c(4,6:10)
ty="wt10"
to="mut10"
slopes[[20]] <- extended3(type,nx,ty,to)
names(slopes)[[20]]<-"dmelanogaster_dilp10"

type = "celegans_mut"

nx=c(1:3,10:12)
ty=9:11
to=6:8
slopes[[21]] <- extended3(type,nx,ty,to,"celegans_mut2")
names(slopes)[[21]]<-"14_wt"

type="dmelanogaster_mut"
nx=c(1:6)
ty="wt10"
to="wt50"
slopes[[22]] <- extended3(type,nx,ty,to)
names(slopes)[[22]]<-"50vs10"


nx=c(7:12)
ty="mut10"
to="mut50"
slopes[[23]] <- extended3(type,nx,ty,to,"dmelanogaster_mut1")
names(slopes)[[23]]<-"50vs10mut"


nx=c(4:6,10:12)
ty="wt50"
to="mut50"
slopes[[24]] <- extended3(type,nx,ty,to,"dmelanogaster_mut2")
names(slopes)[[24]]<-"50mut"


nx=c(1:3,7:9)
ty="wt10"
to="mut10"
slopes[[25]] <- extended3(type,nx,ty,to,"dmelanogaster_mut3")
names(slopes)[[25]]<-"10mut"

#Peforming wilcoxon test for all introns.
slopes<-medn
medn<-sapply(slopes,"[[",1)
medn1<-medn[c(13,12,14,5,6,11,3,4,8)]
medn1<-lapply(medn1,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn1,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn1,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn1,function(x)x$conf.int[[2]]))
ypos <- c(0.1, 0.3, 0.9, 1.3,1.7, 3, 3.4, 3.7, 4)

#Plotting of the results
library(RColorBrewer)
colors <- brewer.pal(9, "Paired")
colors2 <- brewer.pal(3, "Set3")
plot(avg, ypos, xlim=c(-2,2), col= colors[rep(2, 9)], xlab="", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,6), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1.5)

x = ypos
abline(v=0, lty=2,lwd = 1.2)
#text(y = ypos+0.15, x = 1, labels =  c("14 d vs 1 d",  "21 d vs 1 d", "50 d vs 10 d",  "24 mo vs 3.5 mo", "Old vs. young ind","24 mo vs 6 mo",  "Senescent vs Proliferating", "Senescent vs Proliferating","Progeria", "Old vs young ind"), pos = 1, xpd = TRUE, cex = 1.5)
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)

medn2<-medn[c(19,20,15,9,18,10)]
medn2<-lapply(medn2,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn2,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn2,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn2,function(x)x$conf.int[[2]]))
xpos <- c(0.1, 0.9, 1.2, 2, 2.2, 2.4)+0.1
points(avg, xpos, col="#F39200" , xlab="", ylab="",  las=2, pch=19, xlim=c(0,3), lwd = 1.2, cex = 2)
arrows( plus, xpos,  minus, xpos, length=0.03, angle=90, code=3, lwd = 1.2)
medn_wilcox<-medn

#Peforming wilcoxon test for consistent introns.
medn<-sapply(slopes,"[[",2)
names(medn)
medn1<-medn[c(13,12,14,6,11,3,4,8)]
unlist(lapply(medn1,nrow))
medn1<-lapply(medn1,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn1,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn1,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn1,function(x)x$conf.int[[2]]))
ypos <- c(0.1, 0.3, 0.9,1.7, 3, 3.4, 3.7, 4)

library(RColorBrewer)
points(avg, ypos, xlim=c(-1,1), col= colors[rep(2, 9)], xlab="", ylab="", pch=0, ylim=c(0,6), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1.5)

x = ypos
abline(v=0, lty=2,lwd = 1.2)
text(y = ypos+0.15, x = 1, labels =  c("14 d vs 1 d",  "21 d vs 1 d", "50 d vs 10 d",  "24 mo vs 3.5 mo", "Old vs. young ind","24 mo vs 6 mo",  "Senescent vs Proliferating", "Senescent vs Proliferating", "Old vs young ind"), pos = 1, xpd = TRUE, cex = 1.5)
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2,col="red")

medn2<-medn[c(19,20,15,9,16:18,10)]
medn2<-lapply(medn2,function(areao){wilcox.test(areao$SL, conf.int=T)})
avg<-unlist(lapply(medn2,function(x)x$estimate[[1]]))
plus<-unlist(lapply(medn2,function(x)x$conf.int[[1]]))
minus<-unlist(lapply(medn2,function(x)x$conf.int[[2]]))
xpos <- c(0.1, 0.9, 1.2, 1.5, 1.7, 2, 2.2, 2.4)+0.1
points(avg, xpos, col="#F39200" , xlab="", ylab="",  las=2, pch=0, xlim=c(0,3), lwd = 1.2, cex = 2)
text(x = -1, y = xpos+0.15, labels = c("14 d daf2 vs 14 d", "30 d dilp2,3-5 vs 30 d", "50 d dilp 2,3-5 vs 50 d","3.5 mo DR vs 3.5 mo", "5 mo DR vs 5 mo", "16 mo DR vs 16 mo", "26.5 mo DR vs 26.5 mo","26 mo IRS vs 26 mo"), pos = 1, xpd = TRUE, cex = 1.5)
arrows( plus, xpos,  minus, xpos, length=0.03, angle=90, code=3, lwd = 1.2,col="red")

