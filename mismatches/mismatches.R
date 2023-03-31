library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

byregion <- function(x, mut) #Function to combine rnaseqmut output with gtf file
    {
        setkey(x, V1, V2, V3, V4)
        setDT(mut)
        colnames(mut) <- c("chr",   "pos",   "base",   "mut", "mutfw", "mutrev", "mutfwrep", "mutrevrep", "freqfw")
        mut$pos2 <- mut$pos
        byregionintrons <- foverlaps(mut, x, by.x=c("chr", "pos", "pos2"), by.y=c("V1", "V2", "V3"), nomatch=0)
        return(byregionintrons)
    }


genemut <- function(spe, genesf, nb,order=0) #Two files are necessary: rnaseqmut output file and a gtf file with genomic coordinates and gene expression
{   
    files <- Sys.glob(paste(spe, "*.freq", sep = ""), dirmark = FALSE) #listing rnaseqmut files
    files_col<-gsub(spe,"",files) 
    files_col<-gsub(".freq","",files_col) #getting species names
    t <- lapply(files[c(1:length(files))], function(x)fread(x))

    for (x in 1:length(t))
    {
        names(t[[x]]) <- c("chr",   "pos",   "base",   "mut", "mutfw", "mutrev", "mutfwrep", "mutrevrep")
        t[[x]] <- t[[x]][t[[x]]$pos != -1,]
        t[[x]] <- t[[x]][((t[[x]]$mutfw) >= 100) & ((t[[x]]$mutrev) >= 100),]
        setkey(t[[x]], "chr", "pos", "base")
        t[[x]] <- as.data.frame(t[[x]])
        t[[x]]$freqfw <- (t[[x]]$mutfwrep)/(t[[x]]$mutfw+t[[x]]$mutfwrep)
        t[[x]] <- t[[x]][(((t[[x]]$mutfw >= min_reads) & (t[[x]]$mutfwrep == 1) & (t[[x]]$mutrevrep == 0)) | ((t[[x]]$mutrev >= min_reads) & (t[[x]]$mutfwrep == 0) & (t[[x]]$mutrevrep == 1))),]
    }
    for (df in 1:length(t))
        colnames(t[[df]]) <- c("chr",   "pos",   "base", paste("mut", df, sep=""), paste("V5",df,sep=""), paste("V6",df,sep=""), paste("V7",df, sep=""), paste("V8",df, sep=""))
    genesexpr <- fread(paste(spe,"/count_gene_junction_ensembl.bed", sep=""), sep="\t")
    if(length(order)!=1){
      order=c(1:6,6+order)
      genesexpr<-genesexpr[,.SD,.SDcols=order]}
    bygenes <- lapply(t, function(y)byregion(genesexpr, y))
    
    bygenest <- lapply(bygenes, function(x)x[ ,by=c("chr", "V2", "V3", "V4"), .N])
    bygenest %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2, by=c("chr", "V2", "V3", "V4")), .) -> bygenesmutt

    bygenesmutt[is.na(bygenesmutt)] <- 0    
    genesexprbymut <- merge(as.data.frame(genesexpr), bygenesmutt, by.x=c("V1", "V2", "V3"), by.y=c("chr", "V2", "V3"))
    genesexprbymut <- genesexprbymut[apply(genesexprbymut[,c(7:(nb+7-1))], 1, function(x)sum(x > 100) == nb),]
    genesexprbymut <- genesexprbymut[,c((7+nb+1):(7+nb+nb))]/(genesexprbymut[,c(7:(nb+7-1))])
    colnames(genesexprbymut)<-files_col
    return(genesexprbymut)
}

datal <- list("/data/public/cdebes/workspace/scripts/celegans_mut/","/data/public/cdebes/workspace/scripts/celegans/","/data/public/cdebes/workspace/scripts/celegans_polya/","/data/public/cdebes/workspace/scripts/mmusculus_kidney_CR/","/data/public/cdebes/workspace/scripts/mmusculus/","/data/public/cdebes/workspace/scripts/mmusculus_blood/","/data/public/cdebes/workspace/scripts/mmusculus_hypothalamus/", "/data/public/cdebes/workspace/scripts/dmelanogaster_mut/","/data/public/cdebes/workspace/scripts/dmelanogaster_merged/","/data/public/cdebes/workspace/scripts/hsapiens_polya_IMR90/","/data/public/cdebes/workspace/scripts/hsapiens/","/data/public/cdebes/workspace/scripts/hsapiens_prog/","/data/public/cdebes/workspace/scripts/hsapiens_blood/","/data/public/cdebes/workspace/scripts/hsapiens/","/data/public/cdebes/workspace/scripts/hsapiens_IMR90/","/data/public/cdebes/workspace/scripts/hsapiens_HUVECs/","/data/public/cdebes/workspace/scripts/mmusculus_aldr/","/data/public/cdebes/workspace/scripts/hsapiens/",   "/data/public/cdebes/workspace/scripts/rnorvegicus/", "/data/public/cdebes/workspace/scripts/hsapiens_nascent_IMR90/", "/data/public/cdebes/workspace/scripts/hsapiens_nascent/")
datan <- c(15,12,12,8,6,6,7,12,48,4,6,4,12,6,4,4,18,6,6,4,8) #number of samples
ord<-list(c(1,4,2,5,9,8,3,6,7,11,10,12:15),0,0,c(1,3:5,2,6:8),0,c(5,2,1,6,3,4),c(1,2,4,6,7,3,5),c(1,3,4,2,6,5,7:12),0,0,c(1,6,5,3,4,2),0,c(2,1,8,7,3,10,4,9,6,5,11,12),c(1,6,5,3,4,2),0,0,c(1:4,15,14,16:18,5:13),c(1,6,5,3,4,2),0,0,c(6,1,3,8,5,2,7,4)) #order
genes <- c(rep("/data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans.WBcel235.90.gene.bed", 3), rep("/data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.90.gene.bed", 4), rep("/data/public/cdebes/workspace/genomes/dmelanogaster/Drosophila_melanogaster.BDGP6.90.gene.bed", 2), rep("/data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.89.gene.bed", 7), rep("/data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.90.gene.bed", 1),rep("/data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.89.gene.bed", 1),rep("/data/public/cdebes/workspace/genomes/rnorvegicus/Rattus_norvegicus.Rnor_6.0.85.gene.bed", 1),rep("/data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.89.gene.bed", 2))


min_reads = 100
se_list <- mapply(genemut, datal, genes, datan,ord)

celegansday14vsday7polya <-  wilcox.test(rowMeans(se_list[[3]][,c(3,7,11)], na.rm=TRUE),  rowMeans(se_list[[3]][,c(2,6,10)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
celegansday21vsday7polya <-  wilcox.test(rowMeans(se_list[[3]][,c(4,8,12)], na.rm=TRUE),  rowMeans(se_list[[3]][,c(2,6,10)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 
hsapienspolyaIMR90 <-  wilcox.test(rowMeans(se_list[[20]][,c(3,4)], na.rm=TRUE),  rowMeans(se_list[[20]][,c(1,2)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 
celegansday14vsday7 <-  wilcox.test(rowMeans(se_list[[2]][,c(3,7,11)], na.rm=TRUE),  rowMeans(se_list[[2]][,c(2,6,10)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
celegansday21vsday7 <-  wilcox.test(rowMeans(se_list[[2]][,c(4,8,12)], na.rm=TRUE),  rowMeans(se_list[[2]][,c(2,6,10)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
celegansday14vsday1 <-  wilcox.test(rowMeans(se_list[[2]][,c(3,7,11)], na.rm=TRUE),  rowMeans(se_list[[2]][,c(1,5,9)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
celegansday21vsday1 <-  wilcox.test(rowMeans(se_list[[2]][,c(4,8,12)], na.rm=TRUE),  rowMeans(se_list[[2]][,c(1,5,9)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
med_14_daf2 <- wilcox.test(rowMeans(se_list[[1]][,c(7:9)], na.rm=TRUE),  rowMeans(se_list[[1]][,c(1,2,3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
celeganspol2slowvsn21 <- wilcox.test(rowMeans(se_list[[1]][,c(13:15)], na.rm=TRUE),  rowMeans(se_list[[1]][,c(10:12)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
med_14_ama <- wilcox.test(rowMeans(se_list[[1]][,c(4,5,6)], na.rm=TRUE),  rowMeans(se_list[[1]][,c(1,2,3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
medn_hnblood <- wilcox.test(rowMeans(se_list[[13]][,c(3:4,7:8,11:12)], na.rm=TRUE),  rowMeans(se_list[[13]][,c(1:2,5:6,9:10)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #hsoldvsyoungblood
med_mmf <- wilcox.test(rowMeans(se_list[[5]][,c(4:6)], na.rm=TRUE),  rowMeans(se_list[[5]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #mmoldvsyoungkidney
med_rnlf <- wilcox.test(rowMeans(se_list[[19]][,c(4:6)], na.rm=TRUE),  rowMeans(se_list[[19]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #rnoldvsyoungliver
liverimp5f <- wilcox.test(rowMeans(se_list[[17]][,c(4,10,11)], na.rm=TRUE),  rowMeans(se_list[[17]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #mmoldvsyoungkidney5
liverimp16 <- wilcox.test(rowMeans(se_list[[17]][,c(7,13,14)], na.rm=TRUE),  rowMeans(se_list[[17]][,c(12,5,6)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #mmoldvsyoungkidney16
liverimp26 <- wilcox.test(rowMeans(se_list[[17]][,c(18,18)], na.rm=TRUE),  rowMeans(se_list[[17]][,c(15,16,17)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #mmoldvsyoungkidney26
medn_mmblood <- wilcox.test(rowMeans(se_list[[6]][,c(4:6)], na.rm=TRUE),  rowMeans(se_list[[6]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #mmoldvsyoungblood
med_hsf <- wilcox.test(rowMeans(se_list[[11]][,c(2,4,6)], na.rm=TRUE),  rowMeans(se_list[[11]][,c(1,3,5)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)    #hsoldvsyoungnascentIMR90
medn_IMR90 <- wilcox.test(rowMeans(se_list[[15]][,c(1,3)], na.rm=TRUE),  rowMeans(se_list[[15]][,c(2,4)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #hsoldvsyoungtotalRNAIMR90
medn_HUVECs <- wilcox.test(rowMeans(se_list[[16]][,c(1,3)], na.rm=TRUE),  rowMeans(se_list[[16]][,c(2,4)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)#hsapiens_huvecs_nascent
medn_IMR90_nascent <- wilcox.test(rowMeans(se_list[[20]][,c(3,4)], na.rm=TRUE),  rowMeans(se_list[[20]][,c(1,2)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #hsoldvsyoungtotalRNAIMR90
med_hsp <- wilcox.test(rowMeans(se_list[[21]][,c(3,4)], na.rm=TRUE),  rowMeans(se_list[[21]][,c(1,2)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)    #hsoldvsyoungnascentHUVECs
med_hp <- wilcox.test(rowMeans(se_list[[12]][,c(3,4)], na.rm=TRUE),  rowMeans(se_list[[12]][,c(1,2)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)     #hsprogpatientvsind
med_mmcr <- wilcox.test(rowMeans(se_list[[4]][,c(1:4)], na.rm=TRUE),  rowMeans(se_list[[4]][,c(5:8)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)     #mmyoungdrvsyoungkidney
hsoldvsyoungpolyaIMR90 <- wilcox.test(rowMeans(se_list[[3]][,c(4,8,12)], na.rm=TRUE),  rowMeans(se_list[[3]][,c(1,5,9)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 
dmoldvsyoung <- wilcox.test(rowMeans(se_list[[9]][,c(47,48,1)], na.rm=TRUE),  rowMeans(se_list[[9]][,c(13,24,25)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)
brainwt <- wilcox.test(rowMeans(se_list[[8]][,c(4:6)], na.rm=TRUE),  rowMeans(se_list[[8]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 
brainwt2 <- wilcox.test(rowMeans(se_list[[8]][,c(10:12)], na.rm=TRUE),  rowMeans(se_list[[8]][,c(7:9)], na.rm=TRUE), conf.int = TRUE, paired = TRUE)

brainimp10 <- wilcox.test(rowMeans(se_list[[8]][,c(7:9)], na.rm=TRUE),  rowMeans(se_list[[8]][,c(1:3)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 
brainimp50 <- wilcox.test(rowMeans(se_list[[8]][,c(10:12)], na.rm=TRUE),  rowMeans(se_list[[8]][,c(4:6)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) 

brainimp30f <- wilcox.test(rowMeans(se_list[[9]][,c(35:37)], na.rm=TRUE),  rowMeans(se_list[[9]][,c(13,24,25)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #dmdilpyoungvsyoung 
brainimp50f <- wilcox.test(rowMeans(se_list[[9]][,c(11,12,14)], na.rm=TRUE),  rowMeans(se_list[[9]][,c(47,48,1)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #dmdilpoldvsold
med_mmh <- wilcox.test(rowMeans(se_list[[7]][,c(5:7)], na.rm=TRUE),  rowMeans(se_list[[7]][,c(1,2,3,4)], na.rm=TRUE), conf.int = TRUE, paired = TRUE) #dmdilpoldvsold


library(RColorBrewer)
colors <- brewer.pal(4, "Paired")
colors2 <- brewer.pal(3, "Set3")
per_quantile=10
par(mar=c(6,1,1,1), lwd=1.2, cex=1.2)

avg <- c(celegansday14vsday1$estimate[[1]], celegansday21vsday1$estimate[[1]],   brainwt$estimate[[1]], brainwt2$estimate[[1]], med_mmf$estimate[[1]], medn_mmblood$estimate[[1]], med_rnlf$estimate[[1]], med_hsp$estimate[[1]], med_hsf$estimate[[1]], medn_hnblood$estimate[[1]])

plus <- c(celegansday14vsday1$conf.int[[1]], celegansday21vsday1$conf.int[[1]],  brainwt$conf.int[[1]], brainwt2$conf.int[[1]],  med_mmf$conf.int[[1]], medn_mmblood$conf.int[[1]],med_rnlf$conf.int[[1]],  med_hsp$conf.int[[1]], med_hsf$conf.int[[1]], medn_hnblood$conf.int[[1]])

minus <- c(celegansday14vsday1$conf.int[[2]], celegansday21vsday1$conf.int[[2]],   brainwt$conf.int[[2]], brainwt2$conf.int[[2]],  med_mmf$conf.int[[2]], medn_mmblood$conf.int[[2]], med_rnlf$conf.int[[2]], med_hsp$conf.int[[2]], med_hsf$conf.int[[2]], medn_hnblood$conf.int[[2]])
ypos <- c(0.1, 0.5, 1.2, 1.5, 2.4, 3.4, 3.8, 4.6, 5.2, 5.6)                                                                                                                                        
plot(avg, ypos, xlim=c(-0.008,0.01), col= colors[rep(2, 9)], xlab="", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,7), yaxt = 'n', lwd = 1.2, cex = 1.4)

x = ypos
abline(v=0, lty=2,lwd = 1.2)
arrows( plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
mtext("Mismatch level", side = 1, line = 3.2, cex = 1.4)
mtext("Nb of mismatch per gene", side = 1, line = 4.2, cex = 1.4)


avg <- c(med_14_daf2$estimate[[1]], celeganspol2slowvsn21$estimate[[1]], med_14_ama$estimate[[1]], brainimp30f$estimate[[1]], brainimp50f$estimate[[1]], brainimp10$estimate[[1]], brainimp50$estimate[[1]],med_mmcr$estimate[[1]],  liverimp5f$estimate[[1]], liverimp16$estimate[[1]], liverimp26$estimate[[1]])
plus <- c(med_14_daf2$conf.int[[1]],  celeganspol2slowvsn21$conf.int[[1]], med_14_ama$conf.int[[1]], brainimp30f$conf.int[[1]], brainimp50f$conf.int[[1]], brainimp10$conf.int[[1]], brainimp50$conf.int[[1]],med_mmcr$conf.int[[1]], liverimp5f$conf.int[[1]], liverimp16$conf.int[[1]],liverimp26$conf.int[[1]])
minus <- c(med_14_daf2$conf.int[[2]],  celeganspol2slowvsn21$conf.int[[2]], med_14_ama$conf.int[[2]], brainimp30f$conf.int[[2]], brainimp50f$conf.int[[2]], brainimp10$conf.int[[2]], brainimp50$conf.int[[2]],med_mmcr$conf.int[[2]], liverimp5f$conf.int[[2]], liverimp16$conf.int[[2]], liverimp26$conf.int[[2]])
xpos <- c(0.1, 0.4, 0.9, 1.2, 1.5, 2.3, 2.6, 3, 3.4, 3.8, 4.2)
points(avg, xpos, col= "#F39200", xlab="", ylab="",  las=2, pch=19, xlim=c(0,3), lwd = 1.2, cex = 1.4)
arrows( plus, xpos,  minus, xpos, length=0.03, angle=90, code=3, lwd = 1.2)
