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



compute_area <- function(y)
{
    if (is.null(y)){return(NULL)}
    if (ncol(y) == 0){return(NULL)}
    
    area <- sapply(1:ncol(y), function(x){trapz(1:nrow(y), y[,x])})   
    return(area)
}


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



