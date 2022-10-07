### WILCOXON TEST FOR RARE SPLICING ANALYSIS ###

library(ggplot2)
library(RColorBrewer)
library(writexl)
setwd("/cellnet/Isabell/Splicing_old/rare_splicing/out/new_cutoff/")
setwd("/cellnet/Isabell/Splicing_old/rare_splicing/filtered_clusters")
##Load data---------------------------------------------------------------------
files <- list.files(path="/cellnet/Isabell/Splicing_old/rare_splicing/filtered_clusters", pattern="gene_events_0")
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

files2 <- list.files(path="/cellnet/Isabell/Splicing/rare_splicing/out/new_cutoff/", pattern="cluster_events")
nms2 <- gsub(".RDS", "", as.character(files2))
nms2 <- gsub("2022-07-19_cluster_events_", "", as.character(nms2))
nms2 <- gsub("2022-07-21_cluster_events_", "", as.character(nms2))
nms2 <- gsub("2022-07-25_cluster_events_", "", as.character(nms2))
nms2 <- gsub("2022-07-27_cluster_events_", "", as.character(nms2))
nms2 <- gsub("2022-07-31_cluster_events_0.7_", "", as.character(nms2))
files2 <- as.list(files2)
cluster_events <- lapply(files2, FUN=function(x) readRDS(x))
names(cluster_events) <- nms2
cluster_events <- cluster_events[order(names(cluster_events))]

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

filter_samples<-function(x){x<-x[c(2,3,7,8,21,19,27,28,29,30,22,24,23)]}
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
# plotting
png(paste0("/data/public/apapada1/Figures/figure3_new.png"))
plot(avg, ypos, xlim=c(-0.5,0.5), col= cols, xlab="logFC rare splice events", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,14), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
text(avg, ypos, labels=c("14d vs 1d", "21d vs 1d", "14 d ama1 vs 14d", "14d daf2 vs 14d", 
                         "50d vs 10d", "50d RPII215 vs 50d", "24 mo vs. 3.5 mo", "Old vs young ind","26 mo IRS vs 26 mo", "24 mo vs 6 mo",
                         "Senescent vs Proliferating", "Senescent vs Proliferating", "Old vs young ind"),
     cex=1.2, pos=4, offset=4)
dev.off()

pdf(paste0("/data/public/apapada1/Figures/figure3afinal.pdf"))
plot(avg, ypos, xlim=c(-0.5,0.5), col= cols, xlab="logFC rare splice events", ylab="", frame.plot=FALSE, pch=19, ylim=c(0,14), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus, ypos,  minus, ypos, length=0.03, angle=90, code=3, lwd = 1.2)
text(avg, ypos, labels=c("14d vs 1d", "21d vs 1d", "14 d ama1 vs 14d", "14d daf2 vs 14d", 
                         "50d vs 10d", "50d RPII215 vs 50d", "24 mo vs. 3.5 mo", "Old vs young ind","26 mo IRS vs 26 mo", "24 mo vs 6 mo",
                         "Senescent vs Proliferating", "Senescent vs Proliferating", "Old vs young ind"),
     cex=1.2, pos=4, offset=4)
dev.off()
#prepare data for ggplot
wil_gg <- wil_results_gene
wil_gg$group <- c(rep("age",6), rep("lifespan extending treatment",2),rep("age",7), rep("lifespan extending treatment",3), rep("age",8),rep("lifespan extending treatment",2), "age")
wil_gg$color.codes <- ifelse(wil_gg$group == "age", "blue", "orange")
p <- ggplot(data=wil_gg, aes(x=avg, y=seq(nrow(wil_gg)))) +
         geom_point(aes(colour = factor(group))) +
         scale_colour_manual(breaks = unique(as.character(wil_gg$group)), 
                             values = unique(as.character(wil_gg$color.codes))) +
         labs(colour="")+
         xlab("logFC rare splice events") +
         xlim(-0.5,0.5) +
         theme(axis.title.y =element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()) +
         geom_errorbarh(aes(xmin=minus, xmax=plus), 
                         color="grey30", height=0.2) +
         #geom_text(aes(label = rownames(wil_results_gene)),hjust = 0.5,  vjust = -0.7) +
         geom_vline(xintercept=0, linetype="dotted") 

ggsave(plot=p, filename=paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_logFC_gene_0.7_ggplot.png"), height=8, width=8)

### Wilcox Test of logFC, cluster level

wil_clu   <- lapply(cluster_events, FUN = function(x) {
                      apply(as.data.frame(x[,grep("^logFC", colnames(x))]), 2, FUN = function(y){
                       wilcox.test(y[is.finite(y)], conf.int=T)})
             })

# remove comparisons that don't make sense/we don't need
wil_clu$celegans_mut$`logFC.ama1-daf2` <- NULL
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

avg_clu   <- unlist(lapply(wil_clu,function(x) {lapply(x, function(y) y$estimate[[1]])}))
plus_clu  <- unlist(lapply(wil_clu,function(x) {lapply(x, function(y) y$conf.int[[1]])}))
minus_clu <- unlist(lapply(wil_clu,function(x) {lapply(x, function(y) y$conf.int[[2]])}))
pval_clu  <- unlist(lapply(wil_clu,function(x) {lapply(x, function(y) y$p.value[[1]])}))
ypos_clu  <- seq(length(avg_clu))
colors    <- brewer.pal(9, "Paired")

wil_results_clu <- cbind(as.data.frame(avg_clu), as.data.frame(plus_clu), as.data.frame(minus_clu), as.data.frame(pval_clu))

#Plotting
png(paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_logFC_cluster_0.7.png"))
plot(avg_clu, ypos_clu, xlim=c(-1,1), col= colors[rep(2, 9)], xlab="logFC rare splice events", ylab="", frame.plot=FALSE, pch=1, ylim=c(0,29), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1.5)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus_clu, ypos_clu,  minus_clu, ypos_clu, length=0.03, angle=90, code=3, lwd = 1.2)
text(avg, ypos, labels=c("C.elegans 7d vs 1d", "C.elegans 14d vs 1d", "C.elegans 21d vs 1d", "C.elegans 14d vs 7d", "C.elegans 21d vs 7d", "C.elegans 21d vs 14d",
                         "C.elegans ama1 vs WT", "C.elegans daf2 vs WT", "C.elegans polyA 7d vs 1d", "C.elegans polyA 14d vs 1d", "C.elegans polyA 21d vs 1d", "C.elegans polyA 14d vs 7d", "C.elegans polyA 21d vs 7d", "C.elegans polyA 21d vs 14d",
                         "D.melanogaster old vs young", "D.melanogaster RPII old vs young", "D.melanogaster RPII old vs WT old", "D.melanogaster RPII young vs WT young", "D.melanogaster WT old vs young",
                         "HUVECs senescent vs proliferating", "H.sapiens blood old vs young", "IMR90 senescent vs proliferating", "IMR90 polya senescent vs proliferating",
                         "Progeria vs WT", "M.musculus old vs young", "M.musculus blood old vs young", "M.musculus IRS vs WT", "M.musculus DR vs WT", "R.norvegicus old vs young"),
     cex=0.7, pos=4, offset=2)
dev.off()

#prepare data for ggplot
wil_gg <- wil_results_clu
wil_gg$group <- c(rep("age",6), rep("lifespan extending treatment",2),rep("age",7), rep("lifespan extending treatment",3), rep("age",8),rep("lifespan extending treatment",2), "age")
wil_gg$color.codes <- ifelse(wil_gg$group == "age", "blue", "orange")
p1 <- ggplot(data=wil_gg, aes(x=avg_clu, y=seq(nrow(wil_gg)))) +
      geom_point(aes(colour = factor(group))) +
       scale_colour_manual(breaks = unique(as.character(wil_gg$group)), 
                           values = unique(as.character(wil_gg$color.codes))) +
        labs(colour="")+
        xlab("logFC rare splice events") +
        xlim(-1,1) +
        theme(axis.title.y =element_blank(),
               axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
        geom_errorbarh(aes(xmin=minus_clu, xmax=plus_clu), 
                       color="grey30", height=0.2) +
        #geom_text(aes(label = rownames(wil_results_gene)),hjust = 0.5,  vjust = -0.7) +
        geom_vline(xintercept=0, linetype="dotted")

ggsave(plot=p1, filename=paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_logFC_cluster_0.7_ggplot.png"), height=8, width=8)

# save wilcoxon results
wil_results <- list("per gene" = wil_results_gene, "per cluster" = wil_results_clu)
wil_results <- lapply(wil_results, function(x) cbind(" "=rownames(x), x))
write_xlsx(wil_results, paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(),"_wilcox_results_0.7.xlsx"))

##############
# PAIRED WILCOXON TEST BETWEEN CONDITIONS #
# conditions are compared by a paired test between the relative rare events #
##############

### Paired Wilcox Test between conditions, gene level 

wil_pair       <-lapply(gene_events, FUN = function(df) {
  rel.rare <- df[,grep("^rel.rare", colnames(df))]
  iter<-outer(1:ncol(rel.rare), 1:ncol(rel.rare), function(i,j) i)
  iter[!lower.tri(iter,diag=T)] <- 0
  iter<-apply(iter,2,function(x)x[x>0])
  tests<-lapply(iter[-length(iter)],function(a){
    res<-list()
    for(j in 2:length(a)){
      res[[j-1]]<- wilcox.test(rel.rare[,a[1]],rel.rare[,a[j]],paired=T,conf.int=T)
    }
    return(res)
  })
  tests<-unlist(tests,recursive = F)
  nms<-lapply(iter[-length(iter)],function(a){
    nms<-vector()
    for(j in 2:length(a)){
      nms[[j-1]]<- paste(colnames(rel.rare)[a[1]],colnames(rel.rare)[a[j]],sep="-")
    }
    return(nms)
  })
  nms <- unlist(nms, recursive=F)
  names(tests) <- nms
  return(tests)
 }
)

# remove comparisons that don't make sense/we don't need
wil_pair$celegans_mut$`rel.rare.ama1-rel.rare.daf2` <- NULL

# get estimate, conf.int and p-values
avg_pair   <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$estimate[[1]])}))
plus_pair  <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$conf.int[[1]])}))
minus_pair <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$conf.int[[2]])}))
pval_pair  <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$p.value[[1]])}))
ypos_pair  <- seq(length(avg_pair))
colors    <- brewer.pal(9, "Paired")

wil_results_pair <- cbind(as.data.frame(avg_pair), as.data.frame(plus_pair), as.data.frame(minus_pair), as.data.frame(pval_pair))

png(paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_paired_gene_0.07.png"))
plot(avg_clu, ypos_clu, xlim=c(-1,1), col= colors[rep(2, 9)], xlab="Change in rare splice events", ylab="", frame.plot=FALSE, pch=1, ylim=c(0,19), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1.5)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus_clu, ypos_clu,  minus_clu, ypos_clu, length=0.03, angle=90, code=3, lwd = 1.2)
dev.off()

#prepare data for ggplot
wil_gg <- wil_results_clu
wil_gg$group <- c(rep("age",6), rep("lifespan extending treatment",2), rep("age",3),rep("lifespan extending treatment",3), rep("age",2),rep("lifespan extending treatment",2))
wil_gg$color.codes <- ifelse(wil_gg$group == "age", "blue", "orange")
p3 <- ggplot(data=wil_gg, aes(x=avg_clu, y=seq(nrow(wil_gg)))) +
  geom_point(aes(colour = factor(group))) +
  scale_colour_manual(breaks = unique(as.character(wil_gg$group)), 
                      values = unique(as.character(wil_gg$color.codes))) +
  labs(colour="")+
  xlab("Change in rare splice events") +
  xlim(-1,1) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_errorbarh(aes(xmin=minus_clu, xmax=plus_clu), 
                 color="grey30", height=0.2) +
  #geom_text(aes(label = rownames(wil_results_gene)),hjust = 0.5,  vjust = -0.7) +
  geom_vline(xintercept=0, linetype="dotted")

ggsave(plot=p3, filename=paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_paired_gene_0.07_ggplot.png"), height=8, width=8)


### Paired Wilcox Test between conditions, cluster level

wil_pair_clu       <- lapply(cluster_events, FUN = function(df) {
  rel.rare <- df[,grep("^rel.rare", colnames(df))]
  iter<-outer(1:ncol(rel.rare), 1:ncol(rel.rare), function(i,j) i)
  iter[!lower.tri(iter,diag=T)] <- 0
  iter<-apply(iter,2,function(x)x[x>0])
  tests<-lapply(iter[-length(iter)],function(a){
    res<-list()
    for(j in 2:length(a)){
      res[[j-1]]<- wilcox.test(rel.rare[,a[1]],rel.rare[,a[j]],paired=T,conf.int=T)
    }
    return(res)
  })
  tests<-unlist(tests,recursive = F)
  nms<-lapply(iter[-length(iter)],function(a){
    nms<-vector()
    for(j in 2:length(a)){
      nms[[j-1]]<- paste(colnames(rel.rare)[a[1]],colnames(rel.rare)[a[j]],sep="-")
    }
    return(nms)
  })
  nms <- unlist(nms, recursive=F)
  names(tests) <- nms
  return(tests)
}
)

# remove comparisons that don't make sense/we don't need
wil_pair_clu$celegans_mut$`rel.rare.ama1-rel.rare.daf2` <- NULL

# get estimate, conf.int and p-values
avg_pair_clu   <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$estimate[[1]])}))
plus_pair_clu  <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$conf.int[[1]])}))
minus_pair_clu <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$conf.int[[2]])}))
pval_pair_clu  <- unlist(lapply(wil_pair,function(x) {lapply(x, function(y) y$p.value[[1]])}))
ypos_pair_clu  <- seq(length(avg_pair))
colors_clu     <- brewer.pal(9, "Paired")

wil_results_pair <- cbind(as.data.frame(avg_pair_clu), as.data.frame(plus_pair_clu), as.data.frame(minus_pair_clu), as.data.frame(pval_pair_clu))

#Plotting
png(paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_paired_cluster_0.07.png"))
plot(avg_clu, ypos_clu, xlim=c(-1,1), col= colors[rep(2, 9)], xlab="Change in rare splice events", ylab="", frame.plot=FALSE, pch=1, ylim=c(0,19), yaxt = 'n', lwd = 1.2, cex = 2,cex.axis=1.5)
abline(v=0, lty=2,lwd = 1.2)
arrows(plus_clu, ypos_clu,  minus_clu, ypos_clu, length=0.03, angle=90, code=3, lwd = 1.2)
dev.off()

#prepare data for ggplot
wil_gg <- wil_results_clu
wil_gg$group <- c(rep("age",6), rep("lifespan extending treatment",2), rep("age",3),rep("lifespan extending treatment",3), rep("age",2),rep("lifespan extending treatment",2))
wil_gg$color.codes <- ifelse(wil_gg$group == "age", "blue", "orange")
p3 <- ggplot(data=wil_gg, aes(x=avg_clu, y=seq(nrow(wil_gg)))) +
  geom_point(aes(colour = factor(group))) +
  scale_colour_manual(breaks = unique(as.character(wil_gg$group)), 
                      values = unique(as.character(wil_gg$color.codes))) +
  labs(colour="")+
  xlab("Change in rare splice events") +
  xlim(-1,1) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_errorbarh(aes(xmin=minus_clu, xmax=plus_clu), 
                 color="grey30", height=0.2) +
  #geom_text(aes(label = rownames(wil_results_gene)),hjust = 0.5,  vjust = -0.7) +
  geom_vline(xintercept=0, linetype="dotted")

ggsave(plot=p3, filename=paste0("/cellnet/Isabell/Splicing/rare_splicing/wilcox_test/", Sys.Date(), "_wilcox_paired_cluster_0.07_ggplot.png"), height=8, width=8)

