library(plyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

filter <- function(x, rsize=10, nb=4)
{
  return(x[((x$V9 >= rsize) & (x$V7 >= nb)),])
}

spe="/data/public/apapada1/Revisions/h3/trimmed"
spe=getwd()
filepath <- paste(spe, "/*SJ.out.tab", sep="")
cov <- Sys.glob(filepath)
splice <- lapply(1:length(cov), function(x)fread(cov[x], header=F))
splice <- lapply(splice, filter, 10, 4)
splice <- do.call(rbind, splice)
splice <- cbind(splice, n=rep(1, nrow(splice)))
complete <- splice[, lapply(.SD, sum), by = list(V1, V2, V3, V4, V5, V6)]
complete <- complete[complete$n >= length(cov),]

#write.table(complete, file=paste("/data/public/apapada1/4su/SJ.complete1.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(complete, file=paste("/data/public/apapada1/Revisions/h3/SJ.complete.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)

system("cat /data/public/apapada1/Revisions/h3/SJ.complete.bed | awk '{OFS=\"\t\"} {{if ($4 == 1) {print $1, $2, $3, \"+\", $5} else if ($4 == 2) {print $1, $2, $3, \"-\", $5}}}' | bedtools intersect -wa -wb -a stdin -b /data/public/apapada1/Genomes/Homo_sapiens.GRCh38.90.intron.bed > /data/public/apapada1/Revisions/h3/junanot1.bed")

j <- fread("/data/public/apapada1/Revisions/h3/junanot1.bed", header=F,quote="")

j <-  as.data.frame(j)
head(j)
colnames(j) <- c("chrA", "startA", "endA", "strandA", "reads", "chrB", "startB", "endB", "strandB", "id", "score", "strand")
j$length <- j$endB - j$startB
commonStartEnd <- j[(j$strandA == j$strandB) & (j$startA-1 == j$startB) & (j$endA == j$endB),]
#commonStartEnd <- j[(j$strandA == j$strandB) & (j$startA-1 == j$startB) & (j$endA == j$endB-1),]


commonStart <- j[(j$strandA == j$strandB) & (j$startA-1 == j$startB) & (j$endA < j$endB),]

commonEnd <- j[(j$strandA == j$strandB) & (j$startA-1 > j$startB) & (j$endA == j$endB),]

diff <- j[(j$strandA == j$strandB) & (j$startA-1 > j$startB) & (j$endA < j$endB),]

#write.table(commonStartEnd[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/", spe ,"/commonStartEnd.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(commonStart[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/", spe ,"/commonStart.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(commonEnd[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/", spe ,"/commonEnd.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(commonEnd[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/", spe ,"/diff.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)

commonStartEndNoRec <- commonStartEnd[ ! ((commonStartEnd$strandA %in% commonStart$strandA) & (commonStartEnd$startA %in% commonStart$startA)) | ((commonStartEnd$strandA %in% commonEnd$strandA) & (commonStartEnd$endB %in% commonEnd$endB)),]
commonStartEndRec <- rbind(commonStart, commonEnd, diff)


head(commonStartEndNoRec[,c("chrA", "startA", "endA", "strandA",  "reads", "id", "score", "length")])
head(commonStartEndRec[,c("chrA", "startA", "endA", "strandA",  "reads", "id", "score", "length")])


junctions<-rbind.data.frame(commonStartEndRec,commonStartEndNoRec)
write.table(junctions[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file="IntronJunction1.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(junctions[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file="IntronJunction.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

