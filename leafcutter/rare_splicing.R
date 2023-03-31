### WILCOXON TEST FOR RARE SPLICING ANALYSIS ###

library(ggplot2)
library(RColorBrewer)
library(writexl)
setwd("/cellnet/Isabell/Splicing_old/rare_splicing/out/new_cutoff/")
##Load data---------------------------------------------------------------------
files <- list.files(path="/cellnet/Isabell/Splicing_old/rare_splicing/out/new_cutoff/", pattern="gene_events")
nms <- gsub(".RDS", "", as.character(files))
nms <- gsub("2022-07-19_gene_events_", "", as.character(nms))
nms <- gsub("2022-07-21_gene_events_", "", as.character(nms))
nms <- gsub("2022-07-25_gene_events_", "", as.character(nms))
nms <- gsub("2022-07-27_gene_events_", "", as.character(nms))
nms <- gsub("2022-07-31_gene_events_0.7_", "", as.character(nms))
nms<-gsub("2022-08-04_gene_events_0.7_","",as.character(nms))
files <- as.list(files)
gene_events <- lapply(files, FUN=function(x) readRDS(x))
names(gene_events) <- nms
gene_events <- gene_events[order(names(gene_events))]


##############------------------------------------------------------------------
# WILCOXON TEST OF logFC BETWEEN CONDITIONS #
# conditions are compared by testing the logFC #


### Wilcox Test of logFC, gene level

wil  <- lapply(gene_events, FUN = function(x) {
    x<-subset(x,genes!=0)
               apply(as.data.frame(x[,grep("^logFC", colnames(x))]), 2, FUN = function(y){
                 wilcox.test(y[is.finite(y)], conf.int=T)})
             })

# remove comparisons that don't make sense/we don't need
wil$celegans_mut$`logFC.ama1-daf2` <- NULL
#rename
names(wil$dmelanogaster_merged) <-"logFC.Old-Young"
names(wil$hsapiens) <- "logFC.Senescent-Proliferating"
names(wil$hsapiens_blood) <- "logFC.Old-Young"
names(wil$hsapiens_nascent) <- "logFC.Senescent-Proliferating"
names(wil$hsapiens_polya_IMR90) <- "logFC.Senescent-Proliferating"
names(wil$hsapiens_prog) <- "logFC.Progeria-WT"
names(wil$mmusculus) <- "logFC.Old-Young"
names(wil$mmusculus_blood) <- "logFC.Old-Young"
names(wil$mmusculus_hypothalamus) <- "logFC.IRS-WT"
names(wil$mmusculus_kidney_CR) <- "logFC.DR-WT"
names(wil$rnorvegicus) <- "logFC.Old-Young"

filter_samples<-function(x){x<-x[c(2,3,7,8,19,17,25,26,27,29,20,22,21)]}
avg       <- unlist(lapply(wil,function(x) {lapply(x, function(y) y$estimate[[1]])}))
plus      <- unlist(lapply(wil,function(x) {lapply(x, function(y) y$conf.int[[1]])}))
minus     <- unlist(lapply(wil,function(x) {lapply(x, function(y) y$conf.int[[2]])}))
pval      <- unlist(lapply(wil,function(x) {lapply(x, function(y) y$p.value[[1]])}))

colors    <- brewer.pal(9, "Paired")
avg<-filter_samples(avg)
plus<-filter_samples(plus)
minus<-filter_samples(minus)
pval<-filter_samples(pval)
ypos      <- seq(length(avg))
wil_results_gene <- cbind(as.data.frame(avg), as.data.frame(plus), as.data.frame(minus), as.data.frame(pval))
cols<-colors[rep(2,13)]
cols[c(3,4,6,9)]<-"#F39200"

pdf(paste0("/data/public/apapada1/Figures/figure3b.pdf"))
plot(avg, ypos, xlim=c(-0.5,0.5), col= cols, xlab="logFC rare splice events", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,14), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
text(avg, ypos, labels=c("14d vs 1d", "21d vs 1d", "14 d ama1 vs 14d", "14d daf2 vs 14d", 
                         "50d vs 10d", "50d RPII215 vs 50d", "24 mo vs. 3.5 mo", "Old vs young ind","26 mo IRS vs 26 mo", "24 mo vs 6 mo",
                         "Senescent vs Proliferating", "Senescent vs Proliferating", "Old vs young ind"),
     cex=1.2, pos=4, offset=4)
dev.off()
