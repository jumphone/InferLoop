
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




inferloop.calILS<-function(X, Y, r=0, only_pos=FALSE){
    X=X
    Y=Y
    r=r
    only_pos=only_pos
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
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
    ILS[which(is.na(ILS))]=0
    if(only_pos==TRUE){ILS[which(ILS<0)]=0 }
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
        this_v=as.numeric(hash::values(x=h, keys=this_tag))
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
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
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



inferloop.loadSignalNoGap <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   HEADER=HEADER[2:length(HEADER)]
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1), skip=1)$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)), skip=1)
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



inferloop.cicero<-function(indata, used_coords, genome.df, k=50, window=5e+05,sample_num=100){
    indata=as.matrix(indata)
    used_coords=used_coords
    genome.df=genome.df
    window=window
    k=k
    sample_num=sample_num
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
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, k =k, size_factor_normalize = FALSE)
    conns <- .run_cicero_new(cicero_cds, genome.df, window, sample_num)
    return(conns)
    }



######################################
# Archr

archr.computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
  ){
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


archr.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}


determineOverlapCpp <- function(m, overlapCut) {
    .Call('_ArchR_determineOverlapCpp', PACKAGE = 'ArchR', m, overlapCut)
}

rowCorCpp <- function(idxX, idxY, X, Y) {
    .Call('_ArchR_rowCorCpp', PACKAGE = 'ArchR', idxX, idxY, X, Y)
}


inferloop.archrCoA<-function(MAT,VEC,k=50,maxDist=500000,SEED=123){
    MAT=MAT
    VEC=VEC
    SEED=SEED
    set.seed(SEED)
    k = k
    knnIteration = 500
    overlapCutoff = 0.8
    maxDist = maxDist
    scaleTo = 10^4
    log2Norm = TRUE
    threads = 10
    verbose = TRUE
    #################

    rD=VEC
    idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
    knnObj <- archr.computeKNN(data = rD, query = rD[idx,], k = k)
    keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
    knnObj <- knnObj[keepKnn==0,]
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
        rownames(rD)[knnObj[x, ]]
          }) %>% SimpleList
    peak=rownames(MAT)
    bed=inferloop.splitLoop(peak,'-',3)
    bed=as.data.frame(bed)
    colnames(bed)=c('chr','start','end')
    bed$start=as.integer(bed$start)
    bed$end=as.integer(bed$end)
    peakSet <- with(bed, GRanges(chr, IRanges(start, end)))
    peakSet$idx=c(1:nrow(bed))
    peakSummits <- resize(peakSet, 1, "center")
    peakWindows <- resize(peakSummits, 2*maxDist + 1, "center")

    o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
    o <- o[o[,1] != o[,2],]
    o$seqnames <- seqnames(peakSet)[o[,1]]
    o$idx1 <- peakSet$idx[o[,1]]
    o$idx2 <- peakSet$idx[o[,2]]
    o$correlation <- -999.999
    o$Variability1 <- 0.000
    o$Variability2 <- 0.000

    umat=MAT[,unlist(knnObj)]
    gmat=.generate_mean(umat, rep(1:length(knnObj),each=k))
    gS=colSums(gmat)
    groupMat=gmat
    groupMat <- t(t(groupMat) / gS) * scaleTo
        if(log2Norm){
            groupMat <- log2(groupMat + 1)
            }

    CHR=bed[,1]
    chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))

    for(x in seq_along(chri)){
        this_chr=chri[x]
        print(this_chr)
        this_row_index=which(CHR %in% this_chr)
        #Correlations
        idx <- BiocGenerics::which(o$seqnames==chri[x])
        corVals <- rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
        rowVars <- as.numeric(matrixStats::rowVars(groupMat))
        o[idx,]$correlation <- as.numeric(corVals)
        o[idx,]$Variability1 <- rowVars[o[idx,]$idx1]
        o[idx,]$Variability2 <- rowVars[o[idx,]$idx2]
        }
    o$idx1 <- NULL
    o$idx2 <- NULL
    o <- o[!is.na(o$correlation),]
    o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
    o$Pval <- 2*pt(-abs(o$TStat), length(knnObj) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")

    o$VarQuantile1 <- archr.getQuantiles(o$Variability1)
    o$VarQuantile2 <- archr.getQuantiles(o$Variability2)
    mcols(peakSet) <- NULL
    o@metadata$peakSet <- peakSet

    o$chr1=o$seqnames
    o$start1=bed$start[o$queryHits]
    o$end1=bed$end[o$queryHits]
    o$chr2=o$seqnames
    o$start2=bed$start[o$subjectHits]
    o$end2=bed$end[o$subjectHits]
    o$distance = abs((o$start1+o$end1)/2-(o$start2+o$end2)/2)

    return(o)
    }


###################################################

.cart2clock<-function(x,y,circle){
    phi <- (atan(x/y)/2/pi * circle + ifelse(y >= 0, circle,1.5 * circle))
    phi <- phi %% circle
    output=data.frame(rho = sqrt(x * x + y * y), phi = phi)
    output=as.matrix(output)
    return(output)
    }

.clock2cart<-function(rho, phi, circle){
    output=data.frame(x = rho * sin(phi/circle * 2 * pi), y = rho *
        cos(phi/circle * 2 * pi))
    output=as.matrix(output)
    return(output)
    }



inferloop.splitILS<-function(X, Y, r=0){
    X=X
    Y=Y
    r=r
    #######################
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #######################
    X_mean = mean(X)
    Y_mean = mean(Y)
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ######################
    ILS=inferloop.calILS(X,Y,r)
    #######################
    CVEC=as.matrix(.cart2clock(X_delta,Y_delta,360))
    CVEC[which(is.na(CVEC))]=0
    ANGLE=CVEC[,2]
    ##############################
    ILS_FLAG=rep(-1,length(ILS))
    ILS_FLAG[which(ILS>0)]=1
    #############################
    this_order=order(ANGLE)
    ANGLE_order=ANGLE[this_order]
    ILS_FLAG_order=ILS_FLAG[this_order]
    ##############################
    SM_ILS_FLAG_order=smooth.spline(ILS_FLAG_order)$y
    SM_ILS_FLAG_order_fat=c(SM_ILS_FLAG_order[length(SM_ILS_FLAG_order)],
                         SM_ILS_FLAG_order,
                         SM_ILS_FLAG_order[1])
    #################################################
    ILS_FLAG_order_fat=SM_ILS_FLAG_order_fat
    ILS_FLAG_order_fat[which(SM_ILS_FLAG_order_fat>0)]=1
    ILS_FLAG_order_fat[which(SM_ILS_FLAG_order_fat<=0)]=-1
    ######################################################
    FLAG_order_fat=rep(0,length(ILS_FLAG_order_fat))
    i=2
    while(i<length(ILS_FLAG_order_fat)){
        this_before=ILS_FLAG_order_fat[i-1]
        this_after=ILS_FLAG_order_fat[i+1]
        if(this_before!=this_after & FLAG_order_fat[i-1]!=1){FLAG_order_fat[i]=1}
        i=i+1}
    FLAG_order=FLAG_order_fat[2:(length(ILS_FLAG_order)+1)]
    ##############################
    CUT_ANGLE=ANGLE_order[which(FLAG_order>0)]
    ############################
    CLST=rep(1,length(ANGLE))
    if(max(CUT_ANGLE)-min(CUT_ANGLE)>180){
        CLST[which(ANGLE>max(CUT_ANGLE))]=1
        }else{
        CLST[which(ANGLE>max(CUT_ANGLE))]=length(CUT_ANGLE)+1
        }
    i=1
    while(i<length(CUT_ANGLE)){
        this_lw=CUT_ANGLE[i]
        this_up=CUT_ANGLE[i+1]
        CLST[which(ANGLE>this_lw & ANGLE <= this_up)]=i+1
        i=i+1}
    #########################
    OUT=list()
    OUT$ils=ILS
    OUT$clst=CLST
    return(OUT)
    }

#######################################
#20230519


inferloop.calD <-function(X, Y, r=0){
    X=X
    Y=Y
    r=r
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
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
    return(D)
    }


inferloop.inferD<-function(mat, net, r=0,sep='.And.'){
    library(hash)
    r=r # default r=0
    sep=sep
    ###################
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
    out=matrix(0,ncol=1,nrow=nrow(net))
    colnames(out)='D'
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
        z=inferloop.calD(x,y,r=r)
        out[i,]=z
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    print('finished!')
    return(out)
    }





