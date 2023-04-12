setwd("/data/public/apapada1/Revisions/h3")
ny=c("H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2") #sample names

load("gen_diff_slope") #loading slopes
colnames(area_diff)<-c('chr', 'start', 'end', 'strand', 'ensembl_gene_id',"H3_r1","H3_r2","H4_r1","H4_r2","wt_r1","wt_r2") #adding column names
filter_reads <- apply(area[,ny], 1, function(x)((sum(x < 0) == length(ny)))) #filtering out positive slopes
areao <- na.omit(area[filter_reads,])
ty<-c("wt_r1.y","wt_r2.y") #controls
ynf <- rowMeans(areao[,ty])
to<-c("H3_r1.y","H3_r2.y") #H3 overexpression
onf <- rowMeans(areao[,to])
areao$control<-ynf
areao$H3<-onf
to<-c("H4_r1.y","H4_r2.y")
onf <- rowMeans(areao[,to])

#Plot
library(RColorBrewer)
colors <- brewer.pal(4, "Paired")
ypos <- c(0.2, 0.4)
h3<-wilcox.test(log2((1/areao$H3)/(1/areao$control)),conf.int=T)
h4<-wilcox.test(log2((1/areao$H4)/(1/areao$control)),conf.int=T)
minus <- c(h3$conf.int[[1]], h4$conf.int[[1]])
plus <- c(h3$conf.int[[2]], h4$conf.int[[2]])
avg <- c(h3$estimate[[1]], h4$estimate[[1]])

plot(avg, ypos, xlab="", ylab="", yaxt = 'n',ylim=c(0,0.6) ,xlim=c(-0.8,0.2),pch=19, col= c("orange","grey"),frame.plot=FALSE, cex.axis=1.2,cex=2)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
mtext("Log2FC of elongation rate", side = 1, line = 2.0, cex = 1.2)
abline(v=0, lwd=2, lty=2)
text(y = ypos, x = c((avg[1]-0.25),avg[2]), labels =  c( "H3-GFP overexpression","H4-GFP overexpression") ,pos = 1, xpd = TRUE, cex = 1.1)
dev.off()
