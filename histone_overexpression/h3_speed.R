#This script calculates the slopes of the intronic coverage.

library(IRanges)
library(parallel)
library(S4Vectors)
library(XVector)
library(GenomicRanges)
library(rtracklayer)
library(flux)
library(plyr)
library(caTools)

source("/data/public/cdebes/workspace/transcription_aging/premature/utils.R")

size <- read.table("/data/public/cdebes/workspace/genomes/hsapiens/ensembl.genome.size", header=F, sep="")
chrsize <- size[,2]
names(chrsize) <- size[,1]
chrsize<-chrsize

cl <- makeCluster(6, outfile="/data/public/apapada1/Revisions/h3/er.log")

#Rhpc_initialize()
clusterEvalQ(cl, library(IRanges))
clusterEvalQ(cl, library(S4Vectors))
clusterEvalQ(cl, library(XVector))
clusterEvalQ(cl, library(GenomicRanges))
clusterEvalQ(cl, library(rtracklayer))
clusterEvalQ(cl, library(flux))
clusterEvalQ(cl, library(caTools))
clusterExport(cl, ls())

path = "/data/public/apapada1/Revisions/h3/trimmed/"

filepath_fwd <- paste(path, "*Unique.str2.out.bg", sep="")
coverage_fwd <- Sys.glob(filepath_fwd)

filepath_rev <- paste(path, "*Unique.str1.out.bg", sep="")
coverage_rev <- Sys.glob(filepath_rev)
gff_b <- read.bed('/data/public/apapada1/Revisions/h3/IntronJunction.bed')
gff_b<-subset(gff_b,end-start>1000)
cols<-colnames(gff_b)[1:5]
gff_b<-gff_b[,c(1:4,6)]
colnames(gff_b)<-cols

clusterExport(cl, ls())

cov_data_fwd <- extract_coverage(coverage_fwd, chrsize, cl)
cov_data_rev <- extract_coverage(coverage_rev, chrsize, cl)
gen_auc <- lapply(names(chrsize)[1:25], run.chr.slope.strand.spe, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize[1:25], "all")
area_diff <- do.call(rbind, gen_auc)
save(area_diff, file=paste0("/data/public/apapada1/Revisions/h3/gen_diff_slope"))
stopCluster(cl)
