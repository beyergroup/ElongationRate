#Example script for c.elegans

library(rtracklayer)
library(IRanges)
library(parallel)
library(XVector)
library(GenomicRanges)
library(rtracklayer)
library(flux)

source("/data/public/apapada1/Scripts/utils.R")
chrsize <- c(I=15073434,
             II=15379431,
             III=13783801,
             IV=17493839,
             MtDNA=13794,
             V=30934180,
             X=17718943)

cl <- makeCluster(6, outfile="/data/public/apapada1/cluster.log")

clusterEvalQ(cl, library(IRanges))
clusterEvalQ(cl, library(XVector))
clusterEvalQ(cl, library(GenomicRanges))
clusterEvalQ(cl, library(rtracklayer))
clusterEvalQ(cl, library(flux))
clusterExport(cl, ls())

path = "/data/public/apapada1/Revisions/KallistoNew/celegans/"
cedricpath="/data/public/cdebes/workspace/scripts/celegans/"

filepath_fwd <- paste(path, "*Unique.str2.out.bg", sep="")
coverage_fwd <- Sys.glob(filepath_fwd)

filepath_rev <- paste(path, "*Unique.str1.out.bg", sep="")
coverage_rev <- Sys.glob(filepath_rev)

norec<-read.table(paste0(cedricpath,"NoRecIntronJunction.bed"))
rec<-read.table(paste0(cedricpath,"RecIntronJunction.bed"))
gff_b<-rbind.data.frame(rec,norec)
gff_b<-gff_b[,c(1:4,6)]
colnames(gff_b) <- c('chr','start','end','strand', 'id')


clusterExport(cl, ls())
cov_data_fwd <- extract_coverage(coverage_fwd, chrsize, cl)[1:2]
cov_data_rev <- extract_coverage(coverage_rev, chrsize, cl)[1:2]


gen_data<-lapply(names(chrsize),get.cov.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
gen_data<-do.call(rbind,gen_data)
save(gen_data,file=paste0(path,"gen_data.RData"))

gen_auc <- lapply(names(chrsize),run.chr.y.intercept.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
intercept_start<-do.call(rbind,gen_auc)
save(intercept_start,file=paste0(path,"intercept_start.RData"))

gen_auc <- lapply(names(chrsize),run.chr.x.intercept.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
intercept_end<-do.call(rbind,gen_auc)
save(intercept_end,file=paste0(path,"intercept_end.RData"))

gen_auc <- lapply(names(chrsize),run.chr.slope.strand.spe, cov_data_fwd, cov_data_rev, gff_b, chrsize, "all", cl)
slope<-do.call(rbind,gen_auc)
save(slope,file=paste0(path,"slope.RData"))

stopCluster(cl)
