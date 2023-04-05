#Script to find recursively spliced introns

library(plyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "mmusculus"
}
spe=args[1]
s=args[2]

filter <- function(x, rsize=10, nb=4)
{
    return(x[((x$V9 >= rsize) & (x$V7 >= nb)),])
}

filepath <- paste(spe, "/*SJ.out.tab", sep="") #splice junction output of STAR alignment
cov <- Sys.glob(filepath)
splice <- lapply(1:length(cov), function(x)fread(cov[x], header=F))
splice <- lapply(splice, filter, 20, 6)
splice <- do.call(rbind, splice)
splice <- cbind(splice, n=rep(1, nrow(splice)))
complete <- splice[, lapply(.SD, sum), by = list(V1, V2, V3, V4, V5, V6)]
complete <- complete[complete$n >= length(cov),]

write.table(complete, file=paste("/data/public/cdebes/workspace/scripts/", spe ,"/SJ.complete.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)

system(paste("cat /data/public/cdebes/workspace/scripts/",spe,"/SJ.complete.bed | awk '{OFS=\"\t\"} {{if ($4 == 1) {print $1, $2, $3, \"+\", $5} else if ($4 == 2) {print $1, $2, $3, \"-\", $5}}}' | bedtools intersect -wa -wb -a stdin -b /data/public/cdebes/workspace/genomes/",s,"/intron_junction_ensembl.bed > /data/public/cdebes/workspace/scripts/",spe,"/junanot.bed", sep=""))

j <- fread(paste("/data/public/cdebes/workspace/scripts/",spe,"/junanot.bed",sep=""), header=F)
j <-  as.data.frame(j)
head(j)
colnames(j) <- c("chrA", "startA", "endA", "strandA", "reads", "chrB", "startB", "endB", "strandB", "id", "score", "strand")
j$length <- j$endB - j$startB
commonStartEnd <- j[(j$strandA == j$strandB) & (j$startA-1 == j$startB) & (j$endA == j$endB),]

commonStart <- j[(j$strandA == j$strandB) & (j$startA-1 == j$startB) & (j$endA < j$endB),]

commonEnd <- j[(j$strandA == j$strandB) & (j$startA-1 > j$startB) & (j$endA == j$endB),]

diff <- j[(j$strandA == j$strandB) & (j$startA-1 > j$startB) & (j$endA < j$endB),]

commonStartEndNoRec <- commonStartEnd[ ! ((commonStartEnd$strandA %in% commonStart$strandA) & (commonStartEnd$startA %in% commonStart$startA)) | ((commonStartEnd$strandA %in% commonEnd$strandA) & (commonStartEnd$endB %in% commonEnd$endB)),]
commonStartEndRec <- rbind(commonStart, commonEnd, diff)

head(commonStartEndNoRec[,c("chrA", "startA", "endA", "strandA",  "reads", "id", "score", "length")])
head(commonStartEndRec[,c("chrA", "startA", "endA", "strandA",  "reads", "id", "score", "length")])

write.table(commonStartEndNoRec[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/",spe,"/NoRecIntronJunction.bed", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(commonStartEndRec[,c("chrA", "startA", "endA", "strandA",  "id", "score", "strandA")], file=paste("/data/public/cdebes/workspace/scripts/",spe,"/RecIntronJunction.bed", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
