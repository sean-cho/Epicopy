#' Plot segmented files.
#' 
#' \code{epicopy} Plot resulting segments from \code{LRRtoCNA}.
#' 
#' @description Plots segmentation output for individual samples.
#' 
#' @param CNA CNA output from LRRtoCNA
#' @param which_sample Either a numeric/integer to indicate which sample to plot or
#'    a character vector with the sample name.
#' @param threshold Numeric vector indicating where to plot dotted lines to indicate
#'    threshold. Defaults at 0.15.
#' @param ylab Y-axis labels
#' @param xlab X-axis labels
#' @param ylim Y-axis limits
#' @param title Title. Defaults to sample name with Epicopy segmentation string.
#' @param seg_lwd Line width of the segments
#' @param ... Passes arguments to plot to control everything other than segment.
#' 
#' @details Plots the segmentation results from \code{LRRtoCNA}.
#' 
#' @examples
#' 
#' \dontrun{
#' # Load included example data.
#' data(data_vignette)
#' # Plot
#' plot_segments(epi_seg, which_sample = 1)
#' }
#' 
#' @import BiocGenerics Hmisc GenomeInfoDb GenomicRanges IRanges org.Hs.eg.db TxDb.Hsapiens.UCSC.hg19.knownGene
#' @export plot_segments
#' 
#' @seealso \code{\link[Epicopy]{getLRR}} \code{\link[Epicopy]{LRRtoCNA}} 
#'    \code{\link[Epicopy]{export_gistic}} \code{\link[minfi]{RGChannelSet-class}}
#'    \code{\link[minfi]{read.450k.sheet}} \code{\link[minfi]{read.450k.exp}}

plot_segments <- function(CNA, which_sample, threshold = 0.15,
                          ylab = "Copy Number", xlab = 'Genomic Location (kb)',
                          ylim = c(-2,2), title = NULL,
                          seg_lwd = 2, seg_col = c('skyblue3', 'black'),
                          ...){
  
  # Set up chromosome and plot parameters
  chrLengths <- as.numeric(as.data.frame(seqinfo(
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))[1:22, 1])
  chrLengths <- chrLengths/1000
  names(chrLengths) <- rownames(as.data.frame(
    seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)))[1:22]
  
  chrLengthsRun <- cumsum(chrLengths)
  chrLengthsRun <- append(0, chrLengthsRun)
  chromAxisLoc <- rep(0, 22)
  for(i in 1:22){chromAxisLoc[i] = mean(chrLengthsRun[c(i,i+1)])}
  chrAxes <- chrLengthsRun
  chrLengthsRun <- chrLengthsRun[-23]
  names(chrLengthsRun) <- names(chrLengths)
  totalGenome <- sum(chrLengths)
  
  # Set colors
  if(length(seg_col) != 2) stop('Seg_col does not have 2 colors. Only two colors are supported at the moment.\n')
  if(!(.isColor(seg_col))) stop('A non-valid color is found in seg_col.\n')
  colset <- rep(seg_col, 11)
  names(colset) <- names(chrLengths)
  
  # Subset data for samples
  dummy_main <- CNA$output
  colnames(dummy_main) <- c('Sample', 'Chromosome', 'Start', 'End', 
                            'N_Probes', 'Segment_Mean')
  if(!grepl('chr', dummy_main$Chromosome[1])){
    dummy_main$Chromosome <- paste0('chr', dummy_main$Chromosome)
  }
  allsamples <- unique(dummy_main$Sample)
  if(inherits(which_sample, c('numeric', 'integer'))){
    which_sample <- allsamples[which_sample]
  } else
    if(inherits(which_sample, c('character'))){
      if(!any(which_sample %in% allsamples)) 
        stop('Sample name not found in sample list.\n')
    } else {
    stop('Which_sample argument not recognized. See ?plot_segments')
  }
    
  dummy <- subset(dummy_main, Sample %in% which_sample)
  dummy$chrLengths <- chrLengthsRun[as.vector(dummy$Chromosome)]
  dummy$startAdj <- with(dummy, Start/1000 + chrLengths)
  dummy$endAdj <- with(dummy, End/1000 + chrLengths)
  
  if(is.null(title)){
    title <- paste0(unique(dummy$Sample), "\nEpicopy segmentation")
  }
  
  plot(0, xlim = c(0, totalGenome), 
       ylim = ylim,
       type = "n", 
       ylab = ylab,
       xlab = xlab)
  abline(h = 0, col = "grey35", lwd = 1.5)
  abline(v = chrAxes, lty = 2, col = "red")
  axis(side = 3, at = chromAxisLoc, tick = FALSE, labels = c(1:22), cex.axis = 0.5, 
       line= -1)
  with(dummy,
       segments(x0 = startAdj, x1 = endAdj, y0 = Segment_Mean, lwd = seg_lwd,
                col = colset[Chromosome]))
  title(main = title)
  abline(h = c(threshold, -threshold), lwd = 1.5, col = 'grey', lty = 2)
}

.isColor <- function(x)
{
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}

.makecolors <- function (object, pal = cbbPalette, keep.order = TRUE) 
{
  require(RColorBrewer)
  if (is.null(dim(object))) {
    ANS = list()
    CLfactor = inherits(object, "factor")
    if (keep.order & CLfactor) {
      dum.order = levels(object)
    }
    object = as.character(object)
    tempPal = pal
    obj.unique = unique(object)
    obj.colors = tempPal[1:length(obj.unique)]
    names(obj.colors) = obj.unique
    color.set = obj.colors[object]
    names(color.set) = object
    if (keep.order & CLfactor) {
      obj.colors = obj.colors[dum.order]
    }
    ANS$Legend = obj.colors
    ANS$Colors = color.set
    return(ANS)
  }
  if (is.null(dim(object)) == FALSE) {
    color.set.mat = matrix(NA, nrow(object), ncol(object))
    legend.list = list()
    tempPal = pal
    for (i in 1:ncol(object)) {
      object.dum = apply(object, 2, as.character)
      obj.unique = unique(object.dum[, i])
      obj.colors = tempPal[1:length(obj.unique)]
      names(obj.colors) = obj.unique
      color.set.mat[, i] = obj.colors[object.dum[, i]]
      legend.list[[i]] = obj.colors
    }
    colnames(color.set.mat) = colnames(object)
    rownames(color.set.mat) = rownames(object)
    names(legend.list) = colnames(object)
    ANS = list(Legend = legend.list, Colors = color.set.mat)
    return(ANS)
  }
}
