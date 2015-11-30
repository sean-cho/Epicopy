#' @import Biobase minfi

.funnorm = function(rgSet, nPCs = 2, sex = NULL, verbose = TRUE) {
  minfi:::.isRG(rgSet)
  rgSet <- updateObject(rgSet)
  
  if (verbose) 
    cat("[preprocessFunnorm] Mapping to genome\n")
  
  gmSet <- mapToGenome(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)
  
  if (verbose) 
    cat("[preprocessFunnorm] Quantile extraction\n")
  
  extractedData <- minfi:::.extractFromRGSet450k(rgSet)
  
  if (is.null(sex)) {
    gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
    sex <- rep(1L, length(gmSet$predictedSex))
    sex[gmSet$predictedSex == "F"] <- 2L 
  }
  rm(rgSet)
  if (verbose) 
    cat("[preprocessFunnorm] Normalization\n")
  CN <- getCN(gmSet)
  gmSet <- minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData, 
                                         sex = sex, nPCs = nPCs, verbose = subverbose)
  return(gmSet)
  
}
environment(.funnorm) = asNamespace('minfi')


.coerce.pData <- function(pdat) {
  pDat <- apply(pdat, 2, as.character)
  dimnames(pDat) <- dimnames(pdat)
  pDat <- as.data.frame(pDat, stringsAsFactors = FALSE)
  return(pDat)
}