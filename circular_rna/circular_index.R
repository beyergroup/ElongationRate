#Loading the required libraries.

library(dplyr)
library(data.table)

#Here we define the function for the calculation of the circular RNA index.
retrieveCirc <-  function(spe, min_dec)
{
    cat(paste(spe, "\n"))
    files <- Sys.glob(paste("/data/public/cdebes/workspace/scripts/",spe,"/*circ.out.bed", sep = ""), dirmark = FALSE)
    filesSJ <- Sys.glob(paste("/data/public/cdebes/workspace/scripts/",spe,"/*SJ.out.tab_e.bed",sep = ""), dirmark = FALSE)
    t <- lapply(files, read.table)
    tSJ <- lapply(filesSJ, read.table)
    #Count circular RNAs &//
    tgroup <- lapply(t, function(x) x %>% group_by(V1, V2, V3, V6) %>% summarize(CI = sum(V4, na.rm = T)) %>% as.data.frame)
    #Count junctions &//
    tSJgroup1 <- lapply(tSJ, function(x) x %>% group_by(V1, V2, V6) %>% summarize(SJ5 = sum(V4, na.rm = T)) %>% as.data.frame)
    tSJgroup2 <- lapply(tSJ, function(x) x %>% group_by(V1, V3, V6) %>% summarize(SJ3 = sum(V4, na.rm = T)) %>% as.data.frame)
    #Match junction to circular RNA 3' &//
    tcirc <- lapply(1:length(tgroup), function(x)merge(tgroup[[x]], tSJgroup1[[x]], by.x=c("V1", "V3"), by.y=c("V1", "V2")))
    tcirc <- lapply(1:length(tgroup), function(x)merge(tcirc[[x]], tSJgroup2[[x]], by.x=c("V1", "V2"), by.y=c("V1", "V3")))
    #Match junction to circular RNA 5' &//
    tcirc <- lapply(1:length(tgroup), function(x)tcirc[[x]][tcirc[[x]]$CI >= min_dec,])
    for (x in 1:length(tcirc))
    {
        names(tcirc[[x]]) <- c("V1",   "V2",   "V3",   "V6.x", paste("CI",x,sep = ""),   "V6.y", paste("SJ5", x, sep = "_"),  "V6",   paste("SJ3", x, sep = "_"))
    }  
    tcirc %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("V1", "V2", "V3", "V6","V6.x",  "V6.y")), .) -> mergtcirc
    mergtcirc <- mergtcirc[,c(1:3,6,5,7,9:ncol(mergtcirc))]
    tcircratio  <- lapply(seq(5,ncol(mergtcirc),3), function(x)mergtcirc[,x]/(mergtcirc[,x]+((mergtcirc[,x+1]+mergtcirc[,x+2])/2)))
    mergtcirc <- cbind(mergtcirc[,c(1:4)], do.call(cbind, tcircratio))
    return(mergtcirc)
}

list=c("celegans","celegans_mut", "dmelanogaster_mut",  "mmusculus", "mmusculus_blood",  "hsapiens", "hsapiens_nascent", "hsapiens_blood")

tcirc <- lapply(list, function(spe)retrieveCirc(spe, 0))

day14circ <- wilcox.test(rowMeans(tcirc[[1]][,c(7,11,15)]),rowMeans(tcirc[[1]][,c(5,10,13)]), paired = TRUE, conf.int = TRUE)
day14ama1circ <- wilcox.test(rowMeans(tcirc[[2]][,c(4:6)+4]),rowMeans(tcirc[[2]][,c(1:3)+4]), paired = TRUE, conf.int = TRUE)
day21circ <- wilcox.test(rowMeans(tcirc[[1]][,c(8,12,16)]),rowMeans(tcirc[[1]][,c(5,10,13)]), paired = TRUE, conf.int = TRUE)
dmslow10circ <- wilcox.test(rowMeans(tcirc[[3]][,c(7:9)+4]),rowMeans(tcirc[[3]][,c(1:3)+4]), paired = TRUE, conf.int = TRUE)
dmslowcirc <- wilcox.test(rowMeans(tcirc[[3]][,c(10:12)+4]),rowMeans(tcirc[[3]][,c(4:6)+4]), paired = TRUE, conf.int = TRUE)
dmcirc <- wilcox.test(rowMeans(tcirc[[3]][,c(4:6)+4]),rowMeans(tcirc[[3]][,c(1:3)+4]), paired = TRUE, conf.int = TRUE)
mcirc <- wilcox.test(rowMeans(tcirc[[4]][,c(8:10)]),rowMeans(tcirc[[4]][,c(5:7)]), paired = TRUE, conf.int = TRUE)
mbloodcirc <- wilcox.test(rowMeans(tcirc[[5]][,c(8:10)]),rowMeans(tcirc[[5]][,c(5:7)]), paired = TRUE, conf.int = TRUE)
hshuvecstotalcirc <- wilcox.test(rowMeans(tcirc[[6]][,c(6,8,10)]),rowMeans(tcirc[[6]][,c(5,7,9)]), paired = TRUE, conf.int = TRUE)
hsimr90circ <- wilcox.test(rowMeans(tcirc[[7]][,c(7,8)]),rowMeans(tcirc[[7]][,c(5,6)]), paired = TRUE, conf.int = TRUE)
hsbloodcirc <- wilcox.test(rowMeans(tcirc[[8]][,c(7,8,11,12,15,16)]),rowMeans(tcirc[[8]][,c(5,6,9,10,13,14)]), paired = TRUE, conf.int = TRUE)

library(RColorBrewer)
colors <- brewer.pal(4, "Paired")
cols<-colors[rep(2,11)]
cols[c(2,4,5)]<-"#F39200"
                
pdf(paste0("/data/public/apapada1/Figures/Circular.pdf"))
ypos <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4)
minus <- c(day14circ$conf.int[[1]], day14ama1circ$conf.int[[1]],  day21circ$conf.int[[1]], dmslow10circ$conf.int[[1]], dmslowcirc$conf.int[[1]], dmcirc$conf.int[[1]], mcirc$conf.int[[1]], mbloodcirc$conf.int[[1]], hsimr90circ$conf.int[[1]], hshuvecstotalcirc$conf.int[[1]], hsbloodcirc$conf.int[[1]])
plus <- c(day14circ$conf.int[[2]], day14ama1circ$conf.int[[2]],  day21circ$conf.int[[2]], dmslow10circ$conf.int[[2]], dmslowcirc$conf.int[[2]], dmcirc$conf.int[[2]], mcirc$conf.int[[2]], mbloodcirc$conf.int[[2]], hsimr90circ$conf.int[[2]], hshuvecstotalcirc$conf.int[[2]], hsbloodcirc$conf.int[[2]])
avg <- c(day14circ$estimate[[1]], day14ama1circ$estimate[[1]],  day21circ$estimate[[1]], dmslow10circ$estimate[[1]], dmslowcirc$estimate[[1]], dmcirc$estimate[[1]], mcirc$estimate[[1]], mbloodcirc$estimate[[1]], hsimr90circ$estimate[[1]], hshuvecstotalcirc$estimate[[1]], hsbloodcirc$estimate[[1]])
ypos      <- seq(length(avg))
plot(avg*100, ypos, ylab="", yaxt = 'n',xlim=c(-10,18), pch=19, col= cols,frame.plot=FALSE,lwd = 1.2, cex = 2,cex.axis=1, xlab="Circular RNA index")
arrows(plus*100, ypos,  minus*100, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
abline(v=0, lwd=2, lty=2)
text(y = ypos, x = avg*100, labels =  c("14 d vs 1 d", "14 d ama1 vs 14 d", "21 d vs 1 d", "10 d RpII215C4 vs 10 d", "50 d RpII215C4 vs 50 d", "50 d vs 10 d", "24 mo vs 3.5 mo","Old vs young ind", "Senescent vs Proliferating",  "Senescent vs Proliferating", "Old vs Young ind"), cex=1.2, pos=4, offset=4)
dev.off()
