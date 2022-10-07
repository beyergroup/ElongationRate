library(data.table)
cols<-c("chr","start","end","strand","gene_id","score1","chrn","startn","endn","score","fuzziness","overlap")
read_mnase <-function(loc){
    nucleo<-fread(loc, stringsAsFactors = F)
    message("File loaded")
    nucleo<-nucleo[,-5]
    colnames(nucleo)<-cols
    pexon <- nucleo$end - nucleo$start
    nucleo <-  nucleo[pexon > 500,]
    #nucleo <-subset(nucleo,startn!=-1)
    nucleo<-nucleo[order(nucleo$chr,nucleo$startn,nucleo$endn),]
    
    return(nucleo)
}
exon_nucleo_young<-read_mnase("/data/public/cdebes/workspace/scripts/hsapiens_mnase1/exon_nucleo_young.bed")
exon_nucleo_old<-read_mnase("/data/public/cdebes/workspace/scripts/hsapiens_mnase1/exon_nucleo_old.bed")


distFun<-function(nucleoy,nucleoo){
    library(data.table)
     nucleoy[which(nucleoy$overlap<round(nucleoy$endn-nucleoy$startn)/2),12]=0
     nucleoo[which(nucleoo$overlap<round(nucleoo$endn-nucleoo$startn)/2),12]=0
     cat(length(which(nucleoy$overlap==0|(duplicated(nucleoy[,7:9])|duplicated(nucleoy[,7:9],fromLast = T)))))
     nucleoy<-nucleoy[-which(nucleoy$overlap==0&(duplicated(nucleoy[,1:3])|duplicated(nucleoy[,1:3],fromLast = T))),]
     nucleoo<-nucleoo[-which(nucleoo$overlap==0&(duplicated(nucleoo[,1:3])|duplicated(nucleoo[,1:3],fromLast = T))),]
     # 
    setkey(nucleoy, chr, start, end, gene_id)
    nucleo_disty <- nucleoy[,list(dy=mean(diff(startn+((endn-startn)/2))),
                                  c = .N,overlap=mean(overlap)),
                            by=list(chr,start,end,gene_id)]
    
    setkey(nucleoo, chr, start, end,gene_id)
    nucleo_disto <- nucleoo[,list(do=mean(diff(startn+((endn-startn)/2))),
                                  c = .N,overlap=mean(overlap)),
                            by=list(chr,start,end,gene_id)]
    nucleo_dist<-merge(nucleo_disty,nucleo_disto,by=c("chr","start","end","gene_id"))
   
     nucleo_dist[which(nucleo_dist$overlap.x==0),6]<-0
  nucleo_dist[which(nucleo_dist$overlap.y==0),9]<-0
    
        return(nucleo_dist)
}

exons<-distFun(exon_nucleo_young,exon_nucleo_old)
rm(exon_nucleo_old)
rm(exon_nucleo_young)
#save.image("figure4a.RData")
 #load("/data/public/apapada1/R Data/Newer/figure4a.RData")
gc()
intron_nucleo_young<-read_mnase("/data/public/cdebes/workspace/scripts/hsapiens_mnase1/intron_nucleo_young.bed")
intron_nucleo_old<-read_mnase("/data/public/cdebes/workspace/scripts/hsapiens_mnase1/intron_nucleo_old.bed")
gc()
introns<-distFun(intron_nucleo_young,intron_nucleo_old)
rm(intron_nucleo_old)
rm(intron_nucleo_young)
gc()
load("/data/public/apapada1/R Data/Newer/figure4a.RData")
save.image("/data/public/apapada1/Scripts/Figures/Figure4a.RData")
load("/data/public/apapada1/Scripts/Figures/Figure4a.RData")
exons<-subset(exons,end-start<12000)
exons<-subset(exons,c.x>0|c.y>0)
introns<-subset(introns,c.x>0|c.y>0)

exons$densityYoung<-(exons$c.x/(exons$end-exons$start))*1000
exons$densityOld<-(exons$c.y/(exons$end-exons$start))*1000
introns$densityYoung<-(introns$c.x/(introns$end-introns$start))*1000
introns$densityOld<-(introns$c.y/(introns$end-introns$start))*1000

library(RColorBrewer)
colors <- brewer.pal(10, "Paired")



exo<-wilcox.test(exons$densityOld,exons$densityYoung,conf.int=T,paired=T)
int<-wilcox.test(introns$densityOld,introns$densityYoung,conf.int=T,paired=T)
you<-wilcox.test(introns$densityYoung,exons$densityYoung,conf.int=T)
old<-wilcox.test(introns$densityOld,exons$densityOld,conf.int=T)
avg<-c(exo$estimate[[1]],int$estimate[[1]],you$estimate[[1]],old$estimate[[1]])
plus<-c(exo$conf.int[[1]],int$conf.int[[1]],you$conf.int[[1]],old$conf.int[[1]])
minus<-c(exo$conf.int[[2]],int$conf.int[[2]],you$conf.int[[2]],old$conf.int[[2]])
ypos=c(0.1,0.4,0.7,1)
colors <- brewer.pal(4, "Paired")
plot(avg, ypos, col= colors[rep(2, 9)], xlab="", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,1), xlim=c(-0.20,0.02),yaxt = 'n', lwd = 1.2, cex = 1.4)

abline(v=0, lty=2,lwd = 1.2)
text(y = ypos, x = avg, labels=c("Exons\n(senescent vs. proliferating)",  "Introns\n(senescent vs. proliferating)", "Proliferating\n(introns vs. exons)","Senescent\n(introns vs. exons)"), pos = 1, xpd = TRUE, cex = 1)
text(y=0,x=avg[1],labels="NS")
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)

