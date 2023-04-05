run.chr.y.intercept.strand.spe<- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
    #Get annotation
    gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
    gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
    
    #Split by CHR-
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
    
    #Extract genomic region of interest
    sizeofannotf <- nrow(gff_filt_fwd)
    sizeofannotr <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
    {                                                       
        mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
    })
    
    r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
    {                                                       
        mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
    })
    
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
    
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
    
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
    
   cat("Get slope \n")
    gen_slope <- mclapply(gen_data, compute_y_intercept, mc.cores=12)
    
    gen_slope[gen_slope == "NULL"] <- NA
    
    cat(paste(chr, "\n"))
    gen_slope <- do.call(rbind, gen_slope)
    #gen_pvalues <- do.call(rbind, pvalues)
    gen_slope <- cbind(gff_filt, gen_slope)        
    return(gen_slope)
}
get.cov.strand.spe  <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
  #Get annotation
  gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")

  #Split by CHR-
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
  
  #Extract genomic region of interest
  sizeofannotf <- nrow(gff_filt_fwd)
  sizeofannotr <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
  {                                                       
    mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
  })
    
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
  {                                                       
    mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
  })
  
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
  
  gen_data <- c(r_fw, r_rev)

 return(gen_data)
}

run.chr.x.intercept.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl)
{ 
    #Get annotation
    gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
    gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
    
    #if (dim(gff_filt)[1] <= 1)
    #          {return(NULL)}
    
    #Split by CHR-
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
    
    #Extract genomic region of interest
    sizeofannotf <- nrow(gff_filt_fwd)
    sizeofannotr <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
    {                                                       
        mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
    })
    
    r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
    {                                                       
        mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
    })
    
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
    
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
    
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
    
    #Get slope------------------#
    cat("Get slope \n")
    gen_slope <- mclapply(gen_data, compute_x_intercept, mc.cores=12)
    cat(paste(format(object.size(gen_slope), units = "MB"), "\n"))
    
    
    
    cat("alphab")
     gen_slope[gen_slope == "NULL"] <- NA
      cat(paste(chr, "\n"))
    gen_slope <- do.call(rbind, gen_slope)
    gen_slope <- cbind(gff_filt, gen_slope)        
    return(gen_slope)
}


compute_intercept_diff <- function(x1, filt=TRUE)
{            
  filt_val = 0
  read_size = 90
  
  if (is.null(x1))
  {return(NULL)}
  
  slopex <- sapply(1:length(x1), function(z) {sum <- summary(lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]])))); return(sum[4][[1]][[1]])})
  slopey <- sapply(1:length(x1), function(z) {sum <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]]))); return(sum$fit[[1]])})
  slope<-slopex-slopey
  return(slope)                  
}

compute_y_intercept <- function(x1, filt=TRUE)
{
  if (is.null(x1))
  {return(NULL)}
  slope <- sapply(1:length(x1), function(z) {sum <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]]))); return(sum$fit[[1]])})
     return(as.vector(slope))
}

compute_x_intercept <- function(x1, filt=TRUE)
{
  filt_val = 0
  read_size = 90
  
  if (is.null(x1))
  {return(NULL)}

  slope <- sapply(1:length(x1), function(z) {sum <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]]))); return(sum$fit[[length(x1[[z]])]])})
  
 
  return(as.vector(slope))
}

CovIntron <- function(x, n=12)
{
    counts_intron = read.delim(paste("../expression/intron/",x,"/count_intron_junction_ensembl.bed", sep=""), header=F)
    #counts_intron <- counts_intron[apply(counts_intron, 2, function(x)sum(x > 10) > n),]
    counts_intron_minus <- counts_intron[counts_intron$V4 == "-",][,c(1:7,(7+n+1):(7+n*2))]
    counts_intron_minus$intronSize = counts_intron_minus$V3 - counts_intron_minus$V2
    counts_intron_plus <- counts_intron[counts_intron$V4 == "+",][,c(1:(7+n))]
    counts_intron_plus$intronSize = counts_intron_plus$V3 - counts_intron_plus$V2
    colnames(counts_intron_plus)  <-  c('chr', 'start', 'end', 'strand', 'id', 'score',  's', 'y_r1', 'y_r2' , 'y_r3', 'o_r1', 'o_r2', 'o_r3',"intronSize")
    colnames(counts_intron_minus)  <-  c('chr', 'start', 'end', 'strand', 'id', 'score',  's', 'y_r1', 'y_r2' , 'y_r3', 'o_r1', 'o_r2', 'o_r3', "intronSize")
    counts <- rbind(counts_intron_minus, counts_intron_plus)
    cov <- counts[,c(8:13)]/counts$intronSize
    return(cbind(counts[,c(1:7)], cov))
}

CovIntron <- function(x, n=12)
{
    counts_intron = read.delim(paste("../expression/intron/",x,"/count_intron_junction_ensembl.bed", sep=""), header=F)
    #counts_intron <- counts_intron[apply(counts_intron, 2, function(x)sum(x > 10) > n),]
    counts_intron_minus <- counts_intron[counts_intron$V4 == "-",][,c(1:7,20:31)]
    counts_intron_minus$intronSize = counts_intron_minus$V3 - counts_intron_minus$V2
    counts_intron_plus <- counts_intron[counts_intron$V4 == "+",][,c(1:19)]
    counts_intron_plus$intronSize = counts_intron_plus$V3 - counts_intron_plus$V2
    colnames(counts_intron_plus)  <-  c('chr', 'start', 'end', 'strand', 'id', 'score',  's', 'day1_r1', 'day1_r2', 'day1_r3', 'day7_r1', 'day7_r2', 'day7_r3', 'day14_r1', 'day14_r2', 'day14_r3', 'day21_r1', 'day21_r2' , 'day21_r3', "intronSize")
    colnames(counts_intron_minus)  <-  c('chr', 'start', 'end', 'strand', 'id', 'score',  's', 'day1_r1', 'day1_r2', 'day1_r3', 'day7_r1', 'day7_r2', 'day7_r3', 'day14_r1', 'day14_r2', 'day14_r3', 'day21_r1', 'day21_r2' , 'day21_r3', "intronSize")
    counts <- rbind(counts_intron_minus, counts_intron_plus)
    cov <- counts[,c(8:19)]/counts$intronSize
    return(cbind(counts[,c(1:7)], cov))
}



extract_coverage <- function(coverage_filepath, chrSize, cl)
    {
        cov_data <- parLapply(cl, coverage_filepath, function(d)
                 {
                     cat(paste(d, "\n"))
                     chipseq <- import.bedGraph(d)
                     
                     
                     
                     #names of chr cant differ no coverage ...
            
                     
                     seqlengths(chipseq) <- chrSize[levels(seqnames(chipseq))]                       
                     return(chipseq)
                                        #coverage(chipseq,weight=chipseq$score)
                 }
                           )        
        return(cov_data)
    }

#Normfactor <-colSums(sapply(cov_data,function(x)sapply(x,sum,na.rm=T)))
                                        #cov_data <- lapply(1:13,function(x){cov_data[[x]] <- cov_data[[x]]*mean(normFactor)/normFactor[x]}) 

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

getAttributeField <- function (x, field, attrsep = "; ")
{
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = " ", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
        }
        else {
            rv = as.character(NA)
        }
        return(rv)
    })
}

sort_and_remove_overlap <- function(gff)
    {
        #Giving a list of interval from gff files sort and remove overlap // only for non stranded experiments.       
        gff <- gff[((gff$feature == 'exon') | (gff$feature == 'five_prime_UTR') | (gff$feature == 'three_prime_UTR') & (gff$source == 'Coding_transcript')),]
        gff <- gff[order(gff$start),]
        gff <- gff[order(gff$seqname),]
        max <- (dim(gff)[1] - 1)
        gff_overlap <- unique(unlist(mclapply(1:max, function(x) if ((gff[x,]$seqname == gff[x+1,]$seqname) & (gff[x,]$end  >= gff[x+1,]$start) & (gff[x+1,]$end  >= gff[x,]$start)) c(x,x+1), mc.cores=4 )))
        gff_filt <- gff[-c(gff_overlap),]
    }

filter_single_gene <- function(x, gen_elem, direc="")
    {
                                        #Filter for only one exon genes
        if (direc == 'all')
            {
                gff_exon <- x[(x$feature == gen_elem),]
                gene_id <- getAttributeField(gff_exon$attributes, "gene_id")
            }
        
        if (direc == 'minus')
            {
                gff_exon <- x[(x$feature == gen_elem)  & (x$strand == '-'),]
                gene_id <- getAttributeField(gff_exon$attributes, "gene_id")
            }
        
        if (direc == 'plus')
            {
                gff_exon <- x[(x$feature == gen_elem)  & (x$strand == '+'),]
                gene_id <- getAttributeField(gff_exon$attributes, "gene_id")
            }
        
        if (direc == '')
            {
                gff_exon <- x[(x$feature == gen_elem) & (x$source == 'Coding_transcript'),]
                gene_id <- gsub("Parent=Transcript:", "\\1", gff_exon$attributes)
            }

        nb_of_elem_per_gene <- table(gene_id)
        list_single <- names(which(nb_of_elem_per_gene == 1))
        gff_exon$attributes <- gene_id
        list_single <- as.data.frame(list_single)
        dframe_single <- merge(list_single, gff_exon, by.x = "list_single", by.y = "attributes")
        return(dframe_single)
    }

filter_gene <- function(x, gen_elem, chr, chrsizem, direc='', sens='')
    {
        if (direc == 'all')
            {
                gff_exon <- x[(x$chr==chr),]                
            }
        
        
        
        if (sens == 'minus')
            {
                gff_exon <- x[(x$chr==chr) & (x$strand == "-"),]
            }

        if (sens == 'plus')
            {
                gff_exon <- x[(x$chr==chr) & (x$strand == "+"),]
            }

        if (direc == '')
            {
                gff_exon <- x[(x$feature == gen_elem) & (x$source == 'Coding_transcript'),]
            }                          

        return(gff_exon)
    }

extract_genomic_range <- function(x, chr, chrsize, cov_data_chr, cl)
#Extract genomic range of interest
{
    size_cov_data <- length(cov_data_chr)
    l <- dim(x)[1]
    #Rhpc_Export(cl, c('x', 'chrsize'), envir=environment())
    Rhpc_lapply(cl, 1:l, function(g)
                {
                    
                    cat(paste(x[g,]['start'], x[g,]['end'], chr, "\n"))
                    
                    geneCov <- sapply(1:size_cov_data, function(y) 
                                      {                 
                                          chrseq <- cov_data_chr[[y]]
                                          if (is.null(chrseq))
                                              {return(NULL)}    
                                          #Rhpc_Export(cl, c('chrseq'), envir=environment())
                                          
                                        #cat(paste(cov_data_chr[[y]], "\n"))
                                          mapply(function(s,e){window(chrseq, s, e)}, x[g,]$start, x[g,]$end)
                                      }
                                      )
                }
                )
}
    
#    clusterExport(cl, rm(cov_data_chr))
    


extract_genomic_range_irange <- function(x, cov_data_chr, cl)
#Extract genomic range of interest from all the coverage
{
    result <- parSapply(cl, 1:dim(x)[1], function(g)
                     {                         
                                        # small to big%
                         chr <- as.vector(x[g,]$chr)
                                        #chr <- x[g,]$seqname
                         
                         geneCov <- unlist(lapply(1:length(cov_data) ,function(y) 
                                                  {
                                                      chrseq <- cov_data[[y]][[chr]]
                                                      if (is.null(cov_data[[y]][[chr]]))
                                                          {return(NULL)}
                                                      
                                                      unlist(mapply(function(s,e){window(chrseq, s, e)}, x[g,]$start, x[g,]$end), use.names = FALSE)
                                                  }
                                                  )
                                           )
                     }
                        )
    
    return(result)
}




collapse.replicate.rev.minus<- function(x, r, gff_filt, repli)
    {
        #average coverage
        if (is.null(r[[x]]))
            {return(NULL)}
         #if (! length(r[[x]]) == length(r[[x]]))
         #   {return(NULL)}   
        #normFactor<-colSums(a)
        #normFactor[normFactor == 0] <- 1
        #t <- do.call('cbind', lapply(1:12, function(x){a[,x]*mean(normFactor)/normFactor[x]}))
        
        if (gff_filt[x,]$end - gff_filt[x,]$start <= 1)
            {return(NULL)}
        
        if (repli == 1)
            {
                if (gff_filt[x,]$strand == '+')
                    {
                        return(r[[x]])
                    }
                else
                    {
                        return(sapply(r[[x]], rev))
                    }
            }
            
        t <- do.call('cbind', lapply(1:length(r[[x]]), function(y) as.numeric(r[[x]][[y]])))
                        
        if (gff_filt[x,]$strand == "+")
            {
                x1 <- sapply(seq(1,dim(t)[2]-1,repli), function(z) Rle(rowMeans(t[,c(z:(z+repli-1))])))
            }
        else
            {
                x1 <- sapply(seq(1,dim(t)[2]-1,repli), function(z) Rle(rev(rowMeans(t[,c(z:(z+repli-1))]))))
            }
        return(x1)
    }

compute_min <- function(x1, filt=TRUE)
    {                  
        if (is.null(x1))
           {return(NULL)}
        total <- sapply(x1, function(x) sum(x))      
        a <- sapply(1:length(x1), function(z) if (total[z] == 0) {0} else {min(x1[[z]])})
        return(a)                  
    }

compute_avg <- function(x1, filt=TRUE)
    { 
         #if (! length(r[[x]]) == length(r[[x]]))
         #   {return(NULL)}   
        #normFactor<-colSums(a)
        #normFactor[normFactor == 0] <- 1
        #t <- do.call('cbind', lapply(1:12, function(x){a[,x]*mean(normFactor)/normFactor[x]}))
        
        filt_val = 0
        read_size = 90
               
        if (is.null(x1))
            {return(NULL)}
        
        
                                        #Average read control
        
        #total <- sapply(x1, function(x) sum(x)) 
        #avg_read <- sapply(x1, function(x) sum(x)/length(x))
        #if (any(total == 0))
        #    {return(NULL)}
        
        #Absolute value read control over > 9
	#if (is.nan(min(avg_read)))
        #    {return(NULL)}                  
        a <- sapply(1:length(x1), function(z) mean(x1[[z]]))
        return(a)                  
    }

compute_n <- function(x1, filt=TRUE)
{ 
   
    filt_val = 0
    read_size = 90
    
    if (is.null(x1))
    {return(NULL)}
    
    
    a <- sapply(1:200, function(z) mean(x1[[z]]))
    return(a)                  
}



compute_bin_avg <- function(x1, min=10, bin=100, filt=TRUE)
    { 
        total <- sapply(x1, function(x) sum(x)) 
        #avg_read <- sapply(x1, function(x) sum(x)/length(x))
        if (any(total == 0))
            {return(NULL)}
        control_len <- split(as.numeric(x1[[1]]), as.numeric(gl(length(as.numeric(x1[[1]])),round(length(as.numeric(x1[[1]]))/bin),length(as.numeric(x1[[1]])))))
 
        if (length(as.numeric(x1[[1]])) < min)
            {return(NULL)}
        
        a <- sapply(1:length(x1), function(z) sapply(split(as.numeric(x1[[z]]), as.numeric(gl(length(as.numeric(x1[[z]])),round(length(as.numeric(x1[[z]]))/bin),length(as.numeric(x1[[z]]))))), function(y) mean(y)))

        return(a)                  
    }


compute_cf <- function(x1)
    {            
        if (is.null(x1))
           {return(NULL)}
        total <- sapply(x1, function(x) sum(x)) 
        a <- sapply(1:length(x1), function(z) if (total[z] == 0) {as.numeric(cumsum(x1[[z]]))} else {as.numeric(cumsum(x1[[z]])/total[z])}) 
        return(a)                  
    }


compute_slope <- function(x1, filt=TRUE)
    {            
               if (is.null(x1))
           {return(NULL)}
        
        slope <- sapply(1:length(x1), function(z) {sum <- summary(lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]])))); return(sum[4][[1]][[2]])})
        return(slope)                  
    }

compute_slope_sc <- function(x1, filt=TRUE)
{            
        if (is.null(x1))
    {return(NULL)}
    
    slope <- sapply(1:length(x1), function(z) {if(sum(x1[[z]])>0){sum <- summary(lm((scale(as.numeric(x1[[z]]))[,1]) ~ seq(1:length(x1[[z]])))); return(sum[4][[1]][[2]])}else{return(0)}})
    return(slope)                  
}



compute_area <- function(y)
{
    if (is.null(y)){return(NULL)}
    if (ncol(y) == 0){return(NULL)}
    
    area <- sapply(1:ncol(y), function(x){trapz(1:nrow(y), y[,x])})   
    return(area)
}

get_id <- function(x)
    {
        res <- strsplit(as.character(x), "_")[[1]]
        return(paste(res[1], res[2], res[3],sep="_"))
    }

read.bed <- function(file, sepe='\t')
    {
        data <- read.delim(file,header=F, sep=sepe)
        colnames(data) <- c('chr','start','end','strand', 'id')
        data$id <- lapply(data$id, get_id)
        
        #bed <- with(data, GRanges(chr, IRanges(start, end), strand, score, id=id))
        #bed <- renameSeqlevels(bed, c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX', 'chrM')) 
        return(data)
    }

subset.gen.elem <- function(chipseq, gff_b)
    {   
        overlaps <- findOverlaps(gff_b, chipseq)
        idx <- unique(subjectHits(overlaps))
        values <- DataFrame(scores=chipseq$score[idx])
        return(values)
    }

#ranges <- subsetByOverlaps(gff_b, chipseq)
#mcols(ranges) <- c(values)


get.chr.coverage <- function(x, chr, cov_data)
{
    chrs1 <- cov_data[[x]][seqnames(cov_data[[x]]) == chr]
    covs1 <- coverage(chrs1, weight=chrs1$score)[[chr]]
    return(covs1)
}

run.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)               
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                       
                  
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        r <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        cat(paste(format(object.size(r), units = "MB"), "\n"))
        #Avg replicate - reverse minus strand
        #gen_data <- lapply(1:length(r), collapse.replicate.rev.minus, r, gff_filt, repli)
        
	
        ###############################################################################################################
        ## Combine genetic elem / take 1 chr and 10 elem and combine by id give the merge matrix to cumsum function1 ##      
        ###############################################################################################################
        if (opt == "_merged_")
            {
                merge_elem <- function(x)                                                                                 
                    {
                        
                        g <- gen_data[gff_filt$id == x]
                        g[g == "NULL"] <- NA
                        g[is.null(g)] <- NA
                        if (any(is.na(g)))
                            {return(NULL)}
                        
                        f1 <- lapply(1:length(g), function(x)g[[x]][[1]])
                        f2 <- lapply(1:length(g), function(x)g[[x]][[2]])
                        
                        list(do.call(c, f1), do.call(c, f2))
                    }
                
                listid <- unique(unlist(gff_filt$id))
                gen_data <- lapply(listid, merge_elem)
                gff_filt <- data.frame(id=listid)
            }
        
        #out_chr <- list(annot=gff_filt, cov=gen_data)
        #save(out_chr, file=paste("/data/SystemsBiology/cdebes/workspace/transcription_aging/premature/", type, "out_chr_", chr, ref, gen_elem, opt, filtering, mode, ".RData", sep=""))
        #----------------------------#

        
                                        #Get cumsum------------------#
        cat("Get cumsum \n")
        gen_ecdf <- parLapply(cl, r, compute_cf)
        cat("Get AUC \n")       
        gen_auc <- parLapply(cl, gen_ecdf, compute_area)              
        cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

run.chr.cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{
    #Get annotation
    gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
    gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")

    if (dim(gff_filt)[1] <= 1)
              {return(NULL)}
    
    #Split by CHR
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)

    #Extract genomic region of interest
    sizeofannotf <- nrow(gff_filt_fwd)
    sizeofannotr <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
    {                                                       
        mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 32)                                      
    })

    r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
    {                                                       
        mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
    })
    
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
         
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
    
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
    
    #Get cumsum------------------#
    cat("Get cumsum \n")
    gen_ecdf <- mclapply(gen_data, compute_cf, mc.cores=32)
    cat(paste(format(object.size(gen_ecdf), units = "MB"), "\n"))
    
    #Get AUC---------------------#
    cat("Get AUC \n")       
    gen_auc <- mclapply(gen_ecdf, compute_area, mc.cores=32)
    cat(paste(format(object.size(gen_auc), units = "MB"), "\n"))
    cat(paste(chr, "\n"))    
    gen_auc[gen_auc == "NULL"] <- NA
    
    gen_auc <- do.call(rbind, gen_auc)        
    gen_auc <- cbind(gff_filt, gen_auc)        
    return(gen_auc)
}


min.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                       
                  
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        #Avg replicate - reverse minus strand
        #gen_data <- lapply(1:length(r), collapse.replicate.rev.minus, r, gff_filt, repli)
        
	
        ###############################################################################################################
        ## Combine genetic elem / take 1 chr and 10 elem and combine by id give the merge matrix to cumsum function1 ##      
        ########################################################################
    
        #out_chr <- list(annot=gff_filt, cov=gen_data)
        #save(out_chr, file=paste("/data/SystemsBiology/cdebes/workspace/transcription_aging/premature/", type, "out_chr_", chr, ref, gen_elem, opt, filtering, mode, ".RData", sep=""))
        #----------------------------#

        
        #Get cumsum------------------#
        
       
        #cat(gen_ecdf)
        #Get AUC
        gen <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) min(x[y][[1]]))}, mc.cores=32)
        # <- as.vector(gen_ecdf)
        
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        #cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

run.chr.slope.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
    { 
        #Get annotation
        
    gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
    gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
    #Split by CHR
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
    #Extract genomic region of interest
    sizeofannotf <- nrow(gff_filt_fwd)
    sizeofannotr <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
    {                                                       
        mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
    })

    r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
    {                                                       
        mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
    })
    
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
         
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
    
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
    
     #Get slope------------------#
    cat("Get slope \n")
    gen_slope <- mclapply(gen_data, compute_slope, mc.cores=12)
    cat(paste(format(object.size(gen_slope), units = "MB"), "\n"))
        
        
       
        cat("alphab")
        gen_slope[gen_slope == "NULL"] <- NA
       
        cat(paste(chr, "\n"))
        gen_slope <- do.call(rbind, gen_slope)
        gen_slope <- cbind(gff_filt, gen_slope)        
        return(gen_slope)
    }

avg.chr.cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
    #Split by chr
    
    cat(chr)
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
    #Get annotation
    gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
    gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
    cat("alphab")
    
    if (dim(gff_filt)[1] <= 1)
    {return(NULL)}
    
    #Extract genomic region of interest
    sizeofannotf <- nrow(gff_filt_fw)
    sizeofannotr <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
    {                                                       
        mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end)                                      
    })
    
    r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
    {                                                       
        mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 2)       
    })
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
    
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
    
    #Merge vector
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
    #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
    
    
    #Get cumsum------------------#
    
    #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
    gen_auc <- mclapply(gen_data, compute_avg, FALSE,  mc.cores=12)        
    
    # <- as.vector(gen_ecdf)
    cat("alphab")
    #gen_auc <- t(gen_auc)
    #print(data.frame(gff_filt))
    gen_auc[gen_auc == "NULL"] <- NA
    #gen_auc <- do.call(c, gen_auc)
    
    cat(paste(chr, "\n"))
    gen_auc <- do.call(rbind, gen_auc)        
    gen_auc <- cbind(gff_filt, gen_auc)        
    return(gen_auc)
}

avg.chr.n.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
    { 
                                        #Split by chr
       
        cat(chr)
        cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
        cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
        gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
        cat("alphab")
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannotf <- nrow(gff_filt_fw)
        sizeofannotr <- nrow(gff_filt_rev)
        
        r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end)                                      
              })

        r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
              {                                                       
                  mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
               })
        
        if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
         
        if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))

        #Merge vector
        gen_data <- c(r_fw, r_rev)
        gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        gen_auc <- mclapply(gen_data, compute_n, FALSE,  mc.cores=6)        
      
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

avg.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    {        
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                                        
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 10)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        cat(paste(format(object.size(r), units = "MB"), "\n"))
        
        #Get cumsum------------------#
        gen_auc <- parLapply(cl, gen_data, compute_avg, FALSE)
        
      
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }


avg.bin.chr.cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
        cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
        gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
        cat("alphab")
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannotf <- nrow(gff_filt_fw)
        sizeofannotr <- nrow(gff_filt_rev)
        r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end, mc.cores = 32)                                      
              })

        r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
              {                                                       
                  mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
               })
        
        r_fw <- lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))

        rea <- function(x, y)
            {
                
                if (length(x) < 1)
                    {
                        return(NULL)
                    }
                else
                    {   
                        return(x[[y]])
                    }
                
            }
        
        r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, rea, y))
        
        gen_data <- c(r_fw, r_rev)
        gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        gen_auc <- lapply(gen_data, compute_bin_avg, FALSE)
        
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        #gen_auc <- do.call(rbind, gen_auc)        
        #gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

avg.bin.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        
        cat("alphab")
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
              })
                
        gen_data <- r
        
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        gen_auc <- lapply(gen_data, compute_bin_avg, FALSE)
        
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        #gen_auc <- do.call(rbind, gen_auc)        
        #gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }


cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                       
                  
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        cat(paste(format(object.size(r), units = "MB"), "\n"))
        #Avg replicate - reverse minus strand
        #gen_data <- lapply(1:length(r), collapse.replicate.rev.minus, r, gff_filt, repli)
        
	
        ###############################################################################################################
        ## Combine genetic elem / take 1 chr and 10 elem and combine by id give the merge matrix to cumsum function1 ##      
        ########################################################################
    
        #out_chr <- list(annot=gff_filt, cov=gen_data)
        #save(out_chr, file=paste("/data/SystemsBiology/cdebes/workspace/transcription_aging/premature/", type, "out_chr_", chr, ref, gen_elem, opt, filtering, mode, ".RData", sep=""))
        #----------------------------#

        
        #Get cumsum------------------#
        
       
        #cat(gen_ecdf)
        #Get AUC
        gen <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) min(x[y][[1]]))}, mc.cores=32)
        # <- as.vector(gen_ecdf)
        
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        #cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }




avg.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    {        
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                                        
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        cat(paste(format(object.size(r), units = "MB"), "\n"))
        
        #Get cumsum------------------#
        gen_auc <- parLapply(cl, gen_data, compute_avg, FALSE)
             
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }



avg.bin.chr.cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
        cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
        gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
        cat("alphab")
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannotf <- nrow(gff_filt_fw)
        sizeofannotr <- nrow(gff_filt_rev)
        r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end, mc.cores = 32)                                      
              })

        r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
              {                                                       
                  mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
               })
        
        r_fw <- lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))

        rea <- function(x, y)
            {
                
                if (length(x) < 1)
                    {
                        return(NULL)
                    }
                else
                    {   
                        return(x[[y]])
                    }
                
            }
        
        r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, rea, y))
        
        gen_data <- c(r_fw, r_rev)
        gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        gen_auc <- lapply(gen_data, compute_bin_avg, FALSE)
        
        # <- as.vector(gen_ecdf)
        cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        #gen_auc <- do.call(rbind, gen_auc)        
        #gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
     
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                       
                  
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
                
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
           
        return(gen_auc)

    }


max.pos.chr.cov <- function(chr, cov_data, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr <- lapply(1:length(cov_data), get.chr.coverage, chr, cov_data)        
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, direc)
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
        
        
        #Extract genomic region of interest
        sizeofannot <- nrow(gff_filt)
        r <- lapply(1:length(cov_data_chr), function(x)
              {                                       
                  
                  mcmapply(function(s,e){window(cov_data_chr[[x]], s, e)}, gff_filt$start, gff_filt$end, mc.cores = 32)                                      
               })
        
        gen_data <- lapply(1:sizeofannot, function(y)lapply(r, function(x)x[[y]]))
        
        
        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        cat(paste(format(object.size(r), units = "MB"), "\n"))
        #Avg replicate - reverse minus strand
        #gen_data <- lapply(1:length(r), collapse.replicate.rev.minus, r, gff_filt, repli)
        
	
        ###############################################################################################################
        ## Combine genetic elem / take 1 chr and 10 elem and combine by id give the merge matrix to cumsum function1 ##      
        ########################################################################
    
        #out_chr <- list(annot=gff_filt, cov=gen_data)
        #save(out_chr, file=paste("/data/SystemsBiology/cdebes/workspace/transcription_aging/premature/", type, "out_chr_", chr, ref, gen_elem, opt, filtering, mode, ".RData", sep=""))
        #----------------------------#

        
        #Get cumsum------------------#
        
       
        #cat(gen_ecdf)
        #Get AUC
        gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) match(max(x[y][[1]]), x[y][[1]]))}, mc.cores=32)
        # <- as.vector(gen_ecdf)
        
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        #cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
    }

cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)    
        cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
        gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
        cat("alphab")
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannotf <- nrow(gff_filt_fw)
        sizeofannotr <- nrow(gff_filt_rev)
        r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end, mc.cores = 32)                                      
              })

        r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
              {                                                       
                  mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
               })
        
        r_fw <- lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))

        rea <- function(x, y)
            {
                
                if (length(x) < 1)
                    {
                        
                        return(NULL)
                    }
                else
                    {
                        
                        return(x[[y]])
                    }
                
            }
        
        r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, rea, y))
        
        gen_data <- c(r_fw, r_rev)
        gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        #gen_auc <- lapply(gen_data, compute_avg, FALSE)
        
      
        # <- as.vector(gen_ecdf)
        #cat("alphab")
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        #gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        #cat(paste(chr, "\n"))
        #gen_auc <- do.call(rbind, gen_auc)        
        #gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_data)
    }

min.chr.cov.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl)
    { 
        #Split by chr
        cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
        cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
        
        #Get annotation
        gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
        gff_filt_fw <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
        gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
        
        if (dim(gff_filt)[1] <= 1)
            {return(NULL)}
                
        #Extract genomic region of interest
        sizeofannotf <- nrow(gff_filt_fw)
        sizeofannotr <- nrow(gff_filt_rev)
        r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
              {                                                       
                  mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fw$start, gff_filt_fw$end, mc.cores = 32)                                      
              })

        r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
              {                                                       
                  mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 32)       
               })
        
        r_fw <- lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))

        rea <- function(x, y)
            {
                
                if (length(x) < 1)
                    {
                        
                        return(NULL)
                    }
                else
                    {
                        
                        return(x[[y]])
                    }
                
            }
        
        r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, rea, y))
        
        gen_data <- c(r_fw, r_rev)
        gff_filt <- rbind(gff_filt_fw, gff_filt_rev)
                                        #r <- extract_genomic_range(head(gff_filt, 10), chr, chrsize, cov_data_chr, cl)
        

       #Get cumsum------------------#

        #gen_auc <- mclapply(gen_data, function(x){sapply(1:length(x), function(y) (x[y][[1]]))}, mc.cores=32)
        gen_auc <- lapply(gen_data, compute_min, FALSE)
        
      
        # <- as.vector(gen_ecdf)
        
        #gen_auc <- t(gen_auc)
        #print(data.frame(gff_filt))
        gen_auc[gen_auc == "NULL"] <- NA
        #gen_auc <- do.call(c, gen_auc)
        
        cat(paste(chr, "\n"))
        gen_auc <- do.call(rbind, gen_auc)        
        gen_auc <- cbind(gff_filt, gen_auc)        
        return(gen_auc)
}

compute_median <- function(x1, filt=TRUE)
{            
 
  if (is.null(x1))
  {return(NULL)}
  
  #median <- sapply(1:length(x1), function(z) {cumsum <- cumsum(as.numeric(x1[[z]]))/sum(as.numeric(x1[[z]])); median1<- length(cumsum[cumsum<=0.5])/length(cumsum); return(median1)})
 #median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; return(sum(cum@lengths[1:min(which(cum@values>0.5))])/sum(cum@lengths))})
 median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; return(cumsum(cum@lengths)[max(which(cum@values<0.5))]/sum(cum@lengths))})
  #median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; point<-min(which(cum@values>0.5)); return(sum(1:point))})
  return(median)                  
}

compute_median_it <- function(x1, perc, filt=TRUE)
{            
  
  if (is.null(x1))
  {return(NULL)}
  
  median <- sapply(1:length(x1), function(z) {torm<-round(perc*sum(ex@lengths)/100); cumsum <- cumsum(as.numeric(x1[[z]]))/sum(as.numeric(x1[[z]])); median1<- length(cumsum[cumsum<=0.5])/length(cumsum); return(median1)})
  #median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; return(sum(cum@lengths[1:min(which(cum@values>0.5))])/sum(cum@lengths))})
  #median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; return(cumsum(cum@lengths)[max(which(cum@values<0.5))]/sum(cum@lengths))})
  #median <- sapply(1:length(x1), function(z) {cum<-cumsum(x1[[z]])/sum(x1[[z]]) ; point<-min(which(cum@values>0.5)); return(sum(1:point))})
  return(median)                  
}



compute_correlation <- function(x1, filt=TRUE)
{            
  
  if (is.null(x1))
  {return(NULL)}
  
   correlation <- sapply(1:length(x1), function(z) {corr<-cor.test(rep(x1[[z]]@values,x1[[z]]@lengths),1:length(x1[[z]])) ; return(paste(corr[[3]],corr[[4]],sep=":"))})
  return(correlation)                  
}

compute_zeros <- function(x1, filt=TRUE)
{            
  
  if (is.null(x1))
  {return(NULL)}
  
  zeros <- sapply(1:length(x1), function(z) {zeroes<-sum(x1[[z]]@lengths[which(x1[[z]]@values==0)]); return(zeroes/sum(x1[[z]]@lengths))})
  #zeros <- sapply(1:length(x1), function(z) {zeroes<-x1[[z]]; return(zeroes)})
  return(zeros)                  
}



run.chr.median.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
  #Get annotation
  
  gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
  
  #if (dim(gff_filt)[1] <= 1)
  #          {return(NULL)}
  
  #Split by CHR
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
  
  #Extract genomic region of interest
  sizeofannotf <- nrow(gff_filt_fwd)
  sizeofannotr <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
  {                                                       
    mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 2)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
  {                                                       
    mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 2)       
  })
  
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
  
  #Get slope------------------#
  cat("Get slope \n")
  gen_median <- mclapply(gen_data, compute_median, mc.cores=32)
  cat(paste(format(object.size(gen_median), units = "MB"), "\n"))
  
  
  
  cat("alphab")
  #gen_auc <- t(gen_auc)
  #print(data.frame(gff_filt))
  gen_median[gen_median == "NULL"] <- NA
  #slope <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[1,]})
  #pvalues <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[2,]})
  #gen_auc <- do.call(c, gen_auc)        
  cat(paste(chr, "\n"))
  gen_median <- do.call(rbind, gen_median)
  #gen_pvalues <- do.call(rbind, pvalues)
  gen_median <- cbind(gff_filt, gen_median)        
  return(gen_median)
}

run.chr.correlation.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
  #Get annotation
  
  gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
  
  #if (dim(gff_filt)[1] <= 1)
  #          {return(NULL)}
  
  #Split by CHR
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
  
  #Extract genomic region of interest
  sizeofannotf <- nrow(gff_filt_fwd)
  sizeofannotr <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
  {                                                       
    mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 2)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
  {                                                       
    mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 2)       
  })
  
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
  
  #Get slope------------------#
  cat("Get correlation \n")
  gen_correlation <- mclapply(gen_data, compute_correlation, mc.cores=32)
  cat(paste(format(object.size(gen_correlation), units = "MB"), "\n"))
  
  
  
  cat("alphab")
  #gen_auc <- t(gen_auc)
  #print(data.frame(gff_filt))
  gen_correlation[gen_correlation == "NULL"] <- NA
  #slope <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[1,]})
  #pvalues <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[2,]})
  #gen_auc <- do.call(c, gen_auc)        
  cat(paste(chr, "\n"))
  gen_correlation <- do.call(rbind, gen_correlation)
  #gen_pvalues <- do.call(rbind, pvalues)
  gen_correlation <- cbind(gff_filt, gen_correlation)        
  return(gen_correlation)
}



run.chr.zeros.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl)
{ 
  #Get annotation
  
  gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
  
  #if (dim(gff_filt)[1] <= 1)
  #          {return(NULL)}
  
  #Split by CHR
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
  
  #Extract genomic region of interest
  sizeofannotf <- nrow(gff_filt_fwd)
  sizeofannotr <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
  {                                                       
    mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 2)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
  {                                                       
    mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 2)       
  })
  
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
  
  #Get slope------------------#
  cat("Get zeros \n")
  gen_zeros <- mclapply(gen_data, compute_zeros, mc.cores=32)
  cat(paste(format(object.size(gen_zeros), units = "MB"), "\n"))
  
  
  cat("alphab")
  #gen_auc <- t(gen_auc)
  #print(data.frame(gff_filt))
  gen_zeros[gen_zeros == "NULL"] <- NA
  #slope <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[1,]})
  #pvalues <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[2,]})
  #gen_auc <- do.call(c, gen_auc)        
  cat(paste(chr, "\n"))
  gen_zeros <- do.call(rbind, gen_zeros)
  #gen_pvalues <- do.call(rbind, pvalues)
  gen_zeros <- cbind(gff_filt, gen_zeros)        
  return(gen_zeros)
}
# 
# Slidingwindow<-function (data, window, step)
# {
#   total <- length(data)
#   spots <- seq(from = 1, to = (total - window), by = step)
#   result <- list(length = length(spots))
#   for (i in 1:length(spots)) {
#     result[[i]] <- data[spots[i]:(spots[i] +  window - 1)]
#   }
#   return(result)
# }

Slidingwindow<-function (data, window, step) 
   {
       total <- length(data)
       spots <- seq(from = 1, to = total , by = step)
       result <- list(length = length(spots))
       for (i in 1:(length(spots)-1)) {
           result[[i]] <- data[spots[i]:(spots[i] + window - 1)]
         }
       return(result)
     }

compute_sliding <- function(x1, filt=TRUE)
{            
  if (is.null(x1))
  {return(NULL)}
  
 slope <- sapply(1:length(x1), function(z) {sums <- Slidingwindow(as.numeric(x1[[z]]),1000,1000);  vals<-sapply(sums, function(x) {summary(lm(as.numeric(x) ~ seq(1:length(x))))[4][[1]][[2]]}); return(vals)})
 # slope <- sapply(1:length(x1), function(z) {sums <- Slidingwindow(as.numeric(x1[[z]]),1000,1000);  vals<-sapply(sums, function(x) {mean(x)}); return(vals)})
    return(slope)                  
}
 #gc()
# gen_slope <- mclapply(gen_data, compute_sliding, mc.cores=32)


run.chr.sliding.strand.spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, window, cl)
{ 
  #Get annotation
  
  gff_filt <- filter_gene(gff_b, gen_elem, chr, chrsize, "all")
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
  
  #if (dim(gff_filt)[1] <= 1)
  #          {return(NULL)}
  
  #Split by CHR
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get.chr.coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get.chr.coverage, chr, cov_data_rev)
  
  #Extract genomic region of interest
  sizeofannotf <- nrow(gff_filt_fwd)
  sizeofannotr <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x)
  {                                                       
    mcmapply(function(s,e){window(cov_data_chr_fwd[[x]], s, e)}, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 2)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x)
  {                                                       
    mcmapply(function(s,e){rev(window(cov_data_chr_rev[[x]], s, e))}, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 2)       
  })
  
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <-lapply(1:sizeofannotf, function(y)lapply(r_fw, function(x)x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:sizeofannotr, function(y)lapply(r_rev, function(x)x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
  
  #Get slope------------------#
  cat("Get slope \n")
  gen_slope <- mclapply(gen_data, compute_sliding, mc.cores=32)
  gen_slope <- lapply(gen_slope, function(x){as.numeric(x)})
  gen_slope <-lapply(gen_slope, `length<-`, max(lengths(gen_slope)))
  
  cat(paste(format(object.size(gen_slope), units = "MB"), "\n"))
  
  
  
  cat("alphab")
  #gen_auc <- t(gen_auc)
  #print(data.frame(gff_filt))
  gen_slope[gen_slope == "NULL"] <- NA
  #slope <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[1,]})
  #pvalues <- sapply(gen_slope, function(x)if (is.na(x)) {return(NA)} else {x[2,]})
  #gen_auc <- do.call(c, gen_auc)        
  cat(paste(chr, "\n"))
  gen_slope <- do.call(rbind, gen_slope)
  #gen_pvalues <- do.call(rbind, pvalues)
  gen_slope <- cbind(gff_filt, gen_slope)        
  return(gen_slope)
}


