#' Obtain LRR from RGset
#'
#' \code{getLRR} Preprocess RGset and obtain LRR for samples against a 
#'    reference set.
#'  
#' @description Returns a matrix containing the log R ratios (LRR) of samples 
#'    to reference set. Input file has to be an RGset. Normals can be defined 
#'    by the user (see details).
#' 
#' @param rgSet RGChannelSet object containing samples of interest
#' @param Normals Either (1) RGChannelSet of normals (reference set), (2) a 
#'    logical or numeric vector  denoting the normals in the input rgSet, (3) a
#'    character of 'thyroid', 'lung', 'breast', or 'all' denoting reference 
#'    set included in Epicopy, or (4) defaults to NULL, which uses all normals 
#'    in Epicopy.
#' @param sampNames A character vector giving sample names. Defaults to NULL, 
#'    which uses \code{pData(rgSet)$Sample_Name}.
#' @param QN Logical. Perform probe type specific quantile normalization of 
#'    reference set.
#' @param Ref 'mode' or 'median'. Method to calculate reference intensity. 
#'    Defaults to mode using the modeest package.
#' @param mode.bw If Ref = 'mode', uses this bw to estimate the mode. 
#'    Defaults to 0.1.
#' @param mode.method If Ref = 'mode', uses this method to estimate 
#'    the mode. Defaults to 'naive'.
#' @param normal.cnv Logical. Return LRR for reference set.
#' @param mean.center Logical. Perform mean centering of final data.
#' @param filterProbes Logical. Filter probe by TCGA defined 450K probes 
#'    adjacent to a known SNP.
#' @param retainedProbes Character. IF filterProbes is TRUE, denotes probes to keep. 
#'    If NULL, uses probes retained by TCGA.
#' @param keep_fnobj Logical. Save functional normalized object as rda.
#' @param fn_output Directory for fnobject. Defaults to current directory.
#' @param ... Passes argument to preprocessFunnorm.
#' 
#' @details GetLRR takes RGset as an argument as an input.
#' 
#' Normals/Reference samples can be defined multiple ways. Defaults to 
#' NULL which uses all the normals incorporated by the Epicopy package 
#' (see next paragraph). Providing an RGset will BiocGenerics::combine(rgSet, Normals) to 
#' get a dataset. Providing a logical/numeric vector will subset the rgSet 
#' into query and reference, with TRUE/numbers representing reference samples. 
#' Users can also specify one of 'thyroid', 'breast', 'lung', or 'all' to use 
#' Epicopy normals.
#' 
#' Epicopy normals are compiled using solid normal tissue methylation microarrays 
#' from TCGA.
#' 
#' @return Matrix containing LRR values.
#' 
#' @export
#' @import modeest minfi Biobase
#' 
#' @examples
#' 
#' # Get LRR
#' data(data_vignette)
#' # Median is used for quick calculation
#' res_lrr <- getLRR(epi_rg, Normals = 'thyroid', Ref = 'median')
#' 
#' @seealso \code{\link[minfi]{preprocessFunnorm}} 
#'    \code{\link[minfi]{RGChannelSet-class}}
#'

getLRR <- function(rgSet, Normals = NULL, sampNames = NULL, QN = FALSE, 
                   Ref = "mode", mode.bw = 0.1, mode.method = "naive", 
                   normal.cnv = TRUE, mean.center = TRUE, 
                   filterProbes = FALSE, retainedProbes = NULL,
                   keep_fnobj = FALSE, fn_output = NULL, ...) {
  
  # Checks
  
  if (!inherits(rgSet, "RGChannelSet")) {
    stop("rgset has to be RGChannelSet.")
  }
  if (!inherits(Normals, c("numeric", "logical", "RGChannelSet", "NULL", 
                           "character"))) {
    stop("Normals in wrong format.")
  }
  if (!Ref %in% c("mode", "median")) {
    stop("Ref has to be 'mode' or 'median'.")
  }
  
  # Get set for processing depending on Normal input
  
  pData(rgSet) <- .coerce.pData(pData(rgSet))
  
  if(!is.null(Normals) & length(Normals) == 1) if(is.na(Normals)){
    cat('Normals specified as NA. Reference will calculated using all user provided samples.\n')
    Normals <- rep(TRUE, ncol(rgSet))
    normal.cnv <- TRUE
  }
  
  if (inherits(Normals, c("numeric", "logical"))) {
    rgset <- rgSet
    Normals <- Normals
  }
  
  if (inherits(Normals, "RGChannelSet")) {
    if (nrow(rgSet) != nrow(Normals)) 
      stop("Features in normal and tumor set do not match.")
    rgset <- BiocGenerics::combine(rgSet, Normals)
    Normals <- colnames(rgset) %in% colnames(Normals)
  }
  
  if (is.null(Normals)) {
    library(EpicopyData)
    cat("Using all epicopy normals as reference set.\n")
    if (!"all.normals" %in% ls(globalenv())) 
      data("allNormals")
    rgset <- BiocGenerics::combine(rgSet, all.normals)
    Normals <- colnames(rgset) %in% colnames(all.normals)
  }
  
  if (inherits(Normals, "character")) {
    if (!Normals %in% c("all", "thyroid", "breast", "lung")) {
      stop("Character normals have to be one of 'all', 'thyroid', 'breast', or 'lung'.")
    }
    .Normals <- Normals
    if (Normals == "all"){
      library(EpicopyData)
      if(!"all.normals" %in% ls(globalenv())) data('allNormals')
      Normals <- all.normals
    } else
    if (Normals == "thyroid"){
      library(EpicopyData)
      if(!"all.normals" %in% ls(globalenv())) data('allNormals')
      Normals <- all.normals[,pData(all.normals)$Tissue_Type %in% 'thyroid']
    } else
    if (Normals == "breast") {
      library(EpicopyData)
      if(!"all.normals" %in% ls(globalenv())) data('allNormals')
      Normals <- all.normals[,pData(all.normals)$Tissue_Type %in% 'breast']
    } else
    if (Normals == "lung") {
      library(EpicopyData)
      if(!"all.normals" %in% ls(globalenv())) data('allNormals')
      Normals <- all.normals[,pData(all.normals)$Tissue_Type %in% 'lung']
    }
    cat("Using", .Normals, "normals as reference set.\n")
    pData(Normals) <- .coerce.pData(pData(Normals))
    rgset <- BiocGenerics::combine(rgSet, Normals)
    Normals <- colnames(rgset) %in% colnames(Normals)
  }
  
  Normals <- as.logical(Normals)
  if (any(is.na(Normals))) {
    stop("Error in classifying normals. Please review input.")
  }
  
  # Set sample names
  
  if (is.null(sampNames)) {
    colnames(rgset) <- pData(rgset)$Sample_Name
  } else {
    if (ncol(rgset) == length(sampNames)) 
      stop("sampNames is not equal to number of samples in rgset.\n")
    colnames(rgset) <- sampNames
  }
  
  # Funnorm
  
  cat("Performing functional normalization.\n")
  fn.set <- .funnorm(rgset, ...)
  fn.intensity <- getUnmeth(fn.set) + getMeth(fn.set)
  fn.intensity[fn.intensity <= 0] <- 0.01
  fn.intensity <- log2(fn.intensity)
  
  normal.fn <- fn.intensity[, Normals]
  
  if(keep_fnobj) {
    cat('Saving funnorm object with normals to specified directory.')
    if(is.null(fn_output)) fn_output <- '.'
    fn_name <- paste0(fn_output, '/epicopy_fn.rda')
    save(fn.set, file = fn_name)
  }
  
  # Quantile normalization
  
  if (QN) {
    cat("Quantile normalization of normals.\n")
    
    # All probes have been filtered to be autosomal
    normal.t1g <- normal.fn[t1g.probes, ]
    normal.t1r <- normal.fn[t1r.probes, ]
    normal.t2 <- normal.fn[t2.probes, ]
    
    stopifnot(complete.cases(normal.t1g), complete.cases(normal.t1r), 
              complete.cases(normal.t2))
    
    normal.t1g.qn <- normalize.quantiles(normal.t1g)
    dimnames(normal.t1g.qn) <- dimnames(normal.t1g)
    normal.t1r.qn <- normalize.quantiles(normal.t1r)
    dimnames(normal.t1r.qn) <- dimnames(normal.t1r)
    normal.t2.qn <- normalize.quantiles(normal.t2)
    dimnames(normal.t2.qn) <- dimnames(normal.t2)
    normal.fn <- rbind(normal.t1g.qn, normal.t1r.qn, normal.t2.qn)
    
    if (!setequal(autosomal, rownames(normal.fn))) 
      stop("Normalization returns discordant number of probes.")
  }
  
  # Subset into autosomal probe sets
  
  cat("Filtering for autosomal probes.\n")
  normal.fn <- normal.fn[autosomal, ]
  tumor.fn <- fn.intensity[autosomal, !Normals]
  
  # Filter away SNP probes
  
  if (filterProbes) {
    if(is.null(retainedProbes)){
      cat("Filtering for TCGA probeset.")
      normal.fn <- normal.fn[tcga.probeset, ]
      tumor.fn <- tumor.fn[tcga.probeset, ]
    } else {
      if(!inherits(retainedProbes, 'character')) stop('RetainedProbes has to be of
                                                      class character.\n')
      if(!all(retainedProbes %in% rownames(tumor.fn))) stop('Some retainedProbes are
                                                            not found in dataset.\n')
      cat("Filtering for user selected probeset.")
      normal.fn <- normal.fn[retainedProbes,]
      tumor.fn <- tumor.fn[retainedProbes,]
    }
  }
  
  # Calculate reference signal
  
  cat("Calculating reference intensity using", Ref, ".\n")
  if (Ref == "mode") {
    normal.ref <- apply(normal.fn, 1, function(x) {
      mlv(x, bw = mode.bw, method = mode.method, na.rm = TRUE)$M
    })
  }
  if (Ref == "median") {
    normal.ref <- apply(normal.fn, 1, median, na.rm = TRUE)
  }
  
  # Subtract reference
  
  cat("Obtaining LRR.\n")
  tumor.lrr <- sweep(tumor.fn, 1, normal.ref, "-")
  
  if (!normal.cnv) {
    final.lrr <- tumor.lrr
  } else {
    normal.lrr <- sweep(normal.fn, 1, normal.ref, "-")
    final.lrr <- cbind(tumor.lrr, normal.lrr)
  }
  
  # Mean center data
  
  if (mean.center) {
    cat("Mean centering data.\n")
    lrr.mean <- apply(final.lrr, 2, mean)
    final.lrr <- sweep(final.lrr, 2, lrr.mean, "-")
  }
  return(final.lrr)
} 