

library(data.table)

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
            GL=glassoFast::glassoFast(cov_mat, rho_mat,maxIt=1000)

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
    #mean_distance_parameter= 2
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



