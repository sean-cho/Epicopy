#' Convert LRR to segmented copy number
#'
#' \code{LRRtoCNA} Converts LRR to segmented copy number using the circular
#'    binart segmentation (CBS) algorithm.
#'
#' @param LRR Matrix of LRR from \code{getLRR}
#' @param sampID Sample IDs. Default: NULL, uses \code{colnames(LRR)}
#' @param ncores Integer. Default: 1L. Define number of cores for parallel
#'    segmentation. See details.
#' @param seg.undo.splits Undo.splits argument for \code{segment} from
#'    \code{DNAcopy}. Default: \code{'sdundo'}.
#' @param seg.undo.SD undo.SD argument for \code{segment}. Default: 2.
#' @param ... Passes arguments to \code{segment}.
#' 
#' @details This function coerces the LRR file into a suitable input for
#'    the \code{CNA}, \code{smooth.CNA}, and \code{segment} algorithm
#'    provided by the \code{DNAcopy} package. Parallelization for segmentation
#'    of \code{ParDNAcopy} package is implemented when \code{ncores} >= 2L.
#'    
#' @return A segmented DNAcopy object. See \code{\link[DNAcopy]{segment}}.
#' 
#' @import minfi DNAcopy ParDNAcopy
#' @export
#' 
#' @seealso
#'    \code{\link[DNAcopy]{CNA}}
#'    \code{\link[DNAcopy]{smooth.CNA}}
#'    \code{\link[DNAcopy]{segment}}
#'    \code{\link[ParDNAcopy]{parSegment}}
#' 

LRRtoCNA = function(LRR, sampID = NULL, ncores = 1L, 
                    seg.undo.splits = 'sdundo', seg.undo.SD = 2, ...){
  
  if(!all(rownames(LRR) %in% names(hm450))) {stop('Some probe names in LRR does not exist in 450K marker file.')}
  
  # Get marker file (map)
  lrr.map = hm450[rownames(LRR)]
  lrr.chrom = as.numeric(gsub('chr','',seqnames(lrr.map)))
  lrr.loc = lrr.map$probeTarget
  
  # Define sample names
  if(is.null(sampID)){
    sampID = colnames(LRR)
  }
  
  cat('Getting CNA object.\n')
  cna = CNA(genomdat = LRR, chrom = lrr.chrom, maploc = lrr.loc, data.type = 'logratio', sampleid = sampID)
  
  cat('Smoothing CNA object.\n')
  scna = smooth.CNA(cna)
  if(!inherits(ncores,c('numeric','integer'))) {stop('ncores must be numeric/integer')}
  
  ncores = as.integer(ncores)
  if(ncores == 1){
    cat('Running segmentation with 1 core.\n')
    seg = segment(scna, undo.splits = seg.undo.splits, undo.SD = seg.undo.SD, ...)
  } else if(ncores > 1){
    cat('Running segmentation with',ncores,'cores.\n')
    seg = parSegment(CNAobj = scna, undo.splits = seg.undo.splits, undo.SD = seg.undo.SD, njobs = ncores, ...)
  } else {
    stop('ncores is not specified correctly.')
  }
  return(seg)
}