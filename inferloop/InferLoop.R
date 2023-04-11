
library(data.table)
library(stringr)
library(parallel)
library(monocle)
#monocle2, https://www.bioconductor.org/packages/release/bioc/html/monocle.html
library(cicero)
#cicero for monocle2, https://cole-trapnell-lab.github.io/cicero-release/docs/
library(glassoFast)
library(hash)


inferloop.bed2granges<-function(bed){
    library(GenomicRanges)
    colnames(bed)=c('seqname','start','end')
    bed=as.data.frame(bed)
    # split and convert per region
    res <- makeGRangesFromDataFrame(bed)
    return(res)
    }


inferloop.rmOut<-function(X){
    X=X
    Q3=quantile(X,0.75)
    Q1=quantile(X,0.25)
    RANGE=Q3-Q1
    UP=Q3+1.5*RANGE
    LW=Q1-1.5*RANGE
    OUT=X[which(X<UP)]
    OUT=OUT[which(OUT>LW)]
    return(OUT)
    }


inferloop.calABCD<-function(X, Y, X_base, Y_base){
    X=X
    Y=Y
    X_base=X_base
    Y_base=Y_base
    X_delta=X-X_base
    Y_delta=Y-Y_base
    A = sum(X_delta * Y_delta)
    B = sum(X_delta ** 2)
    C = sum(Y_delta ** 2)
    ############################
    D = A / sqrt( B * C )
    OUT=list()
    OUT[['A']]=A
    OUT[['B']]=B
    OUT[['C']]=C
    OUT[['D']]=D
    return(OUT)
    }



inferloop.calILS<-function(X, Y, r=0){
    X=X
    Y=Y
    ############################
    # Ensure positive value
    if(min(X)<0 | min(Y)<0){
        X = X - min(inferloop.rmOut(X))
        Y = Y - min(inferloop.rmOut(Y))
        X[which(X<0)]=0
        Y[which(Y<0)]=0
        }
    ###########################3
    r=r
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = inferloop.calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    D_plus = ( ABCD[['A']] + X_delta * Y_delta ) / sqrt( (ABCD[['B']] + X_delta**2) * (ABCD[['C']] + Y_delta**2) )
    ###########################
    M = D_plus - D
    S = (1-D**2)/(N-1)
    ###########################
    ILS = M / S 
    return( ILS )
    }



inferloop.getUniqLoop <-function(net){
    #######################
    library(hash)
    tmp1=t(as.matrix(net[,c(1,2)]))
    print('sorting ends of each loop...')
    tmp2=apply(tmp1,2,sort)
    tag=apply(tmp2,2,paste0,collapse='.And.')
    utag=unique(tag)
    ####################
    print('hashing...')
    h=hash(keys=utag, values=rep(0,length(utag)))
    ###################
    print('getting unique loop...')
    flag=rep(0,nrow(net))
    i=1
    while(i<=nrow(net)){
        this_tag=tag[i]
        this_v=as.numeric(values(h, this_tag))
        if(this_v>0){flag[i]=1}
        .set(h, keys=this_tag, values=this_v+1)
        if(i %%50000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    out=net[which(flag==0),]
    print(paste0('Number of unique loops: ',nrow(out)))
    print('finished!')
    return(out)
    }



inferloop.inferLoopSignal<-function(mat, net, r=0,sep='.And.'){
    library(hash)
    r=r # default r=0
    sep=sep
    mat=as.matrix(mat)
    net=as.matrix(net)
    tag=rownames(mat)
    net=net[which(net[,1] %in% tag & net[,2] %in% tag),]
    #######################################
    print('hashing...')
    h=hash(keys=tag,values=rep(0,nrow(mat)))
    i=1
    while(i<=nrow(mat)){
        this_tag=tag[i]
        this_v=as.numeric(mat[i,])
        .set(h, keys=this_tag, values=this_v)
        i=i+1}
    #######################################
    print('calculating ILS...')
    out=matrix(0,ncol=ncol(mat),nrow=nrow(net))
    colnames(out)=colnames(mat)
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(values(h, this_tag1)[,1])
        y=as.vector(values(h, this_tag2)[,1])
        z=inferloop.calILS(x,y,r=r)
        out[i,]=z
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    print('finished!')
    return(out)
    }




########################################

inferloop.loadSignal <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1))$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)))
   SIGNAL=as.matrix(SIGNAL)
   rownames(SIGNAL)=NET
   colnames(SIGNAL)=HEADER
   return(SIGNAL)
   }


inferloop.splitLoop <-function(loop,tag='.And.',n=2){
   library(stringr)
   tag=tag
   n=n
   pair=t(matrix(unlist(str_split(loop,tag)),nrow=n))
   return(pair)
   }



inferloop.writeNet <-function(conns, output_path,  cut=400000){
    conns=conns
    CUT=cut
    LOOP=conns[which(rank(-conns[,3], ties.method='random')<= CUT ),]
    LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
    LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
    write.table(LOOP,  file= output_path, row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
    }








#########################################################################################################
.generate_mean <- function(exp_sc_mat, TAG, print_step=50){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG

    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))

    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=rowSums(exp_sc_mat[,this_col])/length(this_col)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }

##############################################################################


inferloop.generateBin <- function(indata, used_coords, n=100, seed=123, type=NULL){
    DATA=indata
    used_coords=used_coords
    if(is.null(type)==FALSE){
        type_coords=as.numeric(as.factor(type))*10
        used_coords=cbind(used_coords, type_coords)
        }
    n=n
    set.seed(seed)
    KM=kmeans(used_coords,centers=n)
    CLST=KM$cluster
    names(CLST)=colnames(indata)
    OUT=.generate_mean(DATA, CLST)
    ##########################
    RETURN=list()
    RETURN$mat=OUT
    RETURN$clst=CLST
    return(RETURN)
    }



inferloop.bin2cell <-function(signal_mat, clst){
    CLST=clst
    SSN=signal_mat
    isLoop=matrix(0,nrow=nrow(SSN),ncol=length(CLST))
    rownames(isLoop)=rownames(SSN)
    colnames(isLoop)=names(CLST)
    UC=unique(CLST)
    i=1
    while(i<=length(UC)){
        this_clst=UC[i]
        this_index1=which(CLST==this_clst)
        this_index2=which(colnames(SSN)==this_clst)
        isLoop[,this_index1]=SSN[,this_index2]
        if(i%%20==1){print(i)}
        i=i+1
        }
    return(isLoop)
    }




######################################
# Cicero

inferloop.getGenomeDF.mm10 <-function(){
    library(GenomicRanges)
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
    genome <- genome[1:21]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }

inferloop.getGenomeDF.hg38 <-function(){
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    genome <- genome[1:24]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }

inferloop.getGenomeDF.hg19 <-function(){
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    genome <- genome[1:24]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }






assemble_connections<-function (cicero_model_list, silent = FALSE){
    types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
    char_hbn <- cicero_model_list[types == "character"]
    gl_only <- cicero_model_list[types == "list"]
    if (!silent) {
        print(paste("Successful cicero models: ", length(gl_only)))
        print("Other models: ")
        print(table(unlist(char_hbn)))
        print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
    }
    cors <- lapply(gl_only, function(gl) {
        cors <- stats::cov2cor(gl$w)
        data.table::melt(as.data.table(cors, keep.rownames = TRUE),
            measure = patterns("[0-9]"))
    })
    cors <- data.table::rbindlist(cors)
    names(cors) <- c("Var1", "Var2", "value")
    data.table::setkey(cors, "Var1", "Var2")
    cors_rec <- as.data.frame(cors[, list(mean_coaccess = mean(value)),
        by = "Var1,Var2"])
    names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
    cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]
    return(cors_rec)
}


generate_windows <- function(window, genomic_coords) {
  if(!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  } else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = window/2)
    l <- r + window - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges=IRanges::IRanges(win_ranges$start,
                                                       win_ranges$end))
  return(gr)
}



get_genomic_range <- function(grs, cds, win) {
  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- cds[(fData(cds)$bp1 < end1 &
                                fData(cds)$bp1 > end2) |
                               (fData(cds)$bp2 < end1 &
                                  fData(cds)$bp2 > end2), ]
  win_range <-
    win_range[as.character(fData(win_range)$chr) ==
                gsub("chr", "",
                     as.character(GenomicRanges::seqnames(grs[win]))),]
  fData(win_range)$mean_bp <-
    (as.numeric(as.character(fData(win_range)$bp1)) +
       as.numeric(as.character(fData(win_range)$bp2)))/2

  return(win_range)
}


calc_dist_matrix <- function(gene_range) {
  dist_mat <- as.matrix(dist(fData(gene_range)$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(fData(gene_range))

  return(dist_mat)
}


get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  return(out)
}



generate_cicero_models<-function (cds, distance_parameter, s = 0.75, window = 5e+05,
    max_elements = 200, genomic_coords = cicero::human.hg19.genome)
{
    assertthat::assert_that(assertthat::is.number(s), s < 1, s > 0)
    grs <- generate_windows(window, genomic_coords)
    fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
    fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
    fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
    print('go')
    outlist <- parallel::mclapply(seq_len(length(grs)), mc.cores = 10,
        function(win) {
            #print(win)
            GL <- "Error"
            win_range <- get_genomic_range(grs, cds, win)
            if (nrow(exprs(win_range)) <= 1) {
                return("Zero or one element in range")
            }
            if (nrow(exprs(win_range)) > max_elements) {
                return("Too many elements in range")
            }
            dist_matrix <- calc_dist_matrix(win_range)
            rho_mat <- get_rho_mat(dist_matrix, distance_parameter,
                s)
            vals <- exprs(win_range)
            cov_mat <- cov(t(vals))
            diag(cov_mat) <- diag(cov_mat) + 1e-04

            #############################################################
            #GL=glasso::glasso(cov_mat, rho_mat)
            GL=glassoFast::glassoFast(cov_mat, rho_mat)

            colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
            colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
            return(GL)
        })
    names_df <- as.data.frame(grs)
    names(outlist) <- paste(names_df$seqnames, names_df$start,
        names_df$end, sep = "_")
    return(outlist)
}


.run_cicero_new<-function(cds,genomic_coords , window=5e+05,sample_num=100){
    cds=cds
    genomic_coords=genomic_coords
    window = window
    silent = FALSE
    sample_num = sample_num
    ########################
    distance_parameters <- estimate_distance_parameter(cds, window = window,
        maxit = 100, sample_num = sample_num, distance_constraint = 250000,
        distance_parameter_convergence = 1e-22, genomic_coords = genomic_coords)
    mean_distance_parameter <- mean(unlist(distance_parameters))
    #####################
    cicero_out <- generate_cicero_models(cds, distance_parameter = mean_distance_parameter,
        window = window, genomic_coords = genomic_coords)
    all_cons <- assemble_connections(cicero_out, silent = silent)
    return(all_cons)
          }



inferloop.cicero<-function(indata, used_coords, genome.df){
    indata=as.matrix(indata)
    used_coords=used_coords
    genome.df=genome.df
    #########################
    library(monocle)
    library(cicero)
    ###########################
    rownames(used_coords)=stringr::str_replace_all(rownames(used_coords),'-','_')
    colnames(indata)=stringr::str_replace_all(colnames(indata),'-','_')
    rownames(indata)=stringr::str_replace_all(rownames(indata),'-','_')
    peakinfo=as.data.frame(t(matrix(unlist(stringr::str_split(rownames(indata),'_')),nrow=3)))
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name=rownames(indata)
    rownames(peakinfo) <- peakinfo$site_name 
    fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
    input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = NULL,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, size_factor_normalize = FALSE)
    conns <- .run_cicero_new(cicero_cds, genome.df)
    return(conns)
    }




















