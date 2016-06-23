#' Plot coverage of 450K across hg19
#' 
#' \code{plotCoverage} Plots the coverage of 450K array across the genome.
#' 
#' @description Plots the coverage of 450K array probes of specified chromosome
#'    range in specified window sizes.
#' 
#' @param chr Chromosome to plot.
#' @param winSize Window size of segments.
#' @param maxVal Upper limit of number of probes for darkest color.
#' @param chrRange Region (in bases) of the chromosome to plot. If NULL plots
#'    the whole chromosome.
#' 
#' @details Uses hg19 annotations. \code{winSize} specifies window sizes and
#'    \code{maxVal} controls the number of probes for the darkest color in the
#'    window.
#' 
#' @import BiocGenerics Hmisc GenomeInfoDb GenomicRanges IRanges
#' @export plotCoverage
#' 
#' @examples
#' 
#' # Plot coverage through the whole of chromosome two
#' plotCoverage('chr2')
#' 
#' # Plot coverage through the whole of chromosome two with windows of 2Mb and max of 200 probes
#' plotCoverage('chr2', winSize = 2e6, maxVal = 200)
#' 

plotCoverage = function(chr, winSize = 1e+06, maxVal = 50, chrRange = NULL, 
                        ...) {
  
  # Set basic vars
  chrLengths = seqlengths(hm450)
  chromosomes = seqnames(hm450)
  
  input = .countProbes(chr = chr, chrRange = chrRange, winSize = winSize)
  if (is.null(chrRange)) {
    from = 1
    to = chrLengths[chr]
  } else {
    if (length(chrRange) != 2) 
      stop("Length of chrRange should be 2")
    chrRange = as.integer(chrRange)
    if (!inherits(chrRange, "integer")) 
      stop("chrRange should be integer")
    if (chrRange[1] < 0) 
      stop("chrRange should be an between 0 and chromosome length")
    if (chrRange[2] > chrLengths[chr]) 
      stop("chrRange should be an between 0 and chromosome length")
    from = chrRange[1]
    to = chrRange[2]
  }
  
  # Calculate plotting parameters in <distances>
  distances = do.call(rbind, lapply(strsplit(names(input), ":", fixed = TRUE), 
                                    unlist))
  distances = as.data.frame(distances)
  distances[, 2:3] = apply(distances[, 2:3], 2, function(x) {
    as.numeric(as.character(x))
  })
  colnames(distances) = c("chr", "start", "end")
  if (is.null(chrRange)) {
    distDiff = (max(chrLengths) - chrLengths[as.character(unique(distances$chr))])/2
  } else {
    distDiff = 0
  }
  
  # Set layout t([1,2,3];[4,4,4])
  layout(matrix(c(1:3, 4, 4, 4), 3, 2), heights = c(0.5, 0.5, 0.25), 
         w = c(0.9, 0.1))
  
  # Plot ideogram
  par(mar = c(0, 2, 3, 2))
  if (is.null(chrRange)) {
    plot(c(0, chrLengths[1]), c(0, 1.8), type = "n", axes = FALSE, xlab = "", 
         ylab = "", xaxs = "i", yaxs = "i", main = paste("Chromosome", 
                                                         gsub("chr", "", unique(distances$chr))))
  } else {
    plot(c(from, to), c(0, 1.8), type = "n", axes = FALSE, xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", main = paste("Chromosome", gsub("chr", 
                                                                 "", unique(distances$chr))))
  }
  print(as.character(unique(distances$chr)))
  bands.ss = subset(bands, chrom == as.character(unique(distances$chr)))
  bands.ss$chromStart = bands.ss$chromStart + distDiff
  bands.ss$chromEnd = bands.ss$chromEnd + distDiff
  apply(bands.ss, 1, .drawIdeogram)
  
  # Plot coverage map
  par(mar = c(4, 2, 0.5, 2))
  if (is.null(chrRange)) {
    plot(0, type = "n", xaxs = "i", yaxs = "i", xlab = NA, ylab = NA, 
         xaxt = "n", yaxt = "n", ylim = c(0, 1), main = NA, axes = FALSE, 
         xlim = c(0, max(chrLengths)))
  } else {
    plot(0, type = "n", xaxs = "i", yaxs = "i", xlab = NA, ylab = NA, 
         xaxt = "n", yaxt = "n", ylim = c(0, 1), main = NA, axes = FALSE, 
         xlim = c(from, to))
  }
  distCol = .makeSeqCol(input, 10, pal = c("white", "blue"), maxVal = maxVal)
  distances = cbind(distances, Color = distCol$colors)
  distances$start = distances$start + distDiff
  distances$end = distances$end + distDiff
  apply(distances, 1, function(x) {
    rect(x[2], 0, x[3], 1, lwd = 0.5, col = x[4])
  })
  axis(side = 1, at = distances[round(seq(1, length(input), length.out = 5)), 2], 
       labels = round(distances[round(seq(1, length(input), length.out = 5)), 2]/1e+06))
  title(xlab = "Genomic Distances (MB)", line = 2)
  rect(0 + distDiff, 0, chrLengths[as.character(unique(distances$chr))] + 
         distDiff, 1, lwd = 1.5)
  
  # Plot legend for coverage map
  par(mar = c(1, 0, 1, 0))
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ylab = NA, xlab = NA)
  legend.image = as.raster(matrix(distCol$legend, nrow = 1))
  rasterImage(legend.image, 0.25, 0.25, 0.75, 1)
  rect(0.25, 0.25, 0.75, 1)
  text(x = 0.25, y = 0.575, pos = 2, "0")
  text(x = 0.75, y = 0.575, pos = 4, paste(">", maxVal, sep = ""))
  winLabel = format(signif(distances$start[2] - distances$start[1], digits = 3), 
                    trim = TRUE, scientific = TRUE)
  winLabel = strsplit(winLabel, "e+", fixed = TRUE)[[1]]
  winMetPrefix = switch(winLabel[2], 
                        `04` = paste(as.character(as.numeric(winLabel[1]) * 10), "kb", sep = ""), 
                        `05` = paste(as.character(as.numeric(winLabel[1]) * 100), "kb", sep = ""), 
                        `06` = paste(as.character(as.numeric(winLabel[1]) * 1), "Mb", sep = ""), 
                        `07` = paste(as.character(as.numeric(winLabel[1]) * 10), "Mb", sep = ""), 
                        `08` = paste(as.character(as.numeric(winLabel[1]) * 100), "Mb", sep = ""))
  mtext(paste("Probe density (#probes/", winMetPrefix, ")", sep = ""), 
        3, line = 0, cex = 0.7)
  axis(side = 1, at = seq(0.25, 0.75, length.out = 5)[2:4], labels = NA, 
       cex.axis = 0.5, line = -0.25, lwd = 0, lwd.ticks = 1, tcl = -0.25)
  axis(side = 1, at = seq(0.25, 0.75, length.out = 5)[2:4], 
       labels = seq(0, maxVal, length.out = 5)[2:4], cex.axis = 0.5, line = -1.2, lwd = 0)
  rect(0.8, 0.25, 0.805, 1, xpd = TRUE)
  text(0.805, 0.625, labels = paste("= ", winMetPrefix, sep = ""), pos = 4)
  
  # Plot legend for ideogram
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ylab = NA, xlab = NA)
  legend("center", title = "Ideogram legend", cex = 0.75, legend = c("Centromere", 
                                                                     "Giemsa Positive 25", "Giemsa Positive 50", "Giemsa Positive 75", 
                                                                     "Giemsa Positive 100", "Giemsa Negative", "Giemsa Variable", "Stalk"), 
         fill = c("maroon", "grey75", "grey60", "grey40", "grey15", "white", 
                  "skyblue2", "tan"))
}

.countProbes = function(chr, chrRange = NULL, winSize = 1e+06, ...) {
  
  # Set basic vars
  chrLengths = seqlengths(hm450)
  chromosomes = seqnames(hm450)
  
  # Checks
  winSize = as.integer(winSize)
  if (!inherits(winSize, "integer")) 
    stop("Window size has to be an integer between 1kb to 10Mb")
  if (winSize >= 1e+07) 
    stop("Window size cannot be larger than 10Mb")
  if (winSize <= 10000) 
    stop("Window size cannot be smaller 1kb")
  if (is.null(chrRange)) {
    from = 1
    to = chrLengths[chr]
  } else {
    if (length(chrRange) != 2) 
      stop("Length of chrRange should be 2")
    chrRange = as.integer(chrRange)
    if (!inherits(chrRange, "integer")) 
      stop("chrRange should be integer")
    if (chrRange[1] < 0) 
      stop("chrRange should be an between 0 and chromosome length")
    if (chrRange[2] > chrLengths[chr]) 
      stop("chrRange should be an between 0 and chromosome length")
    from = chrRange[1]
    to = chrRange[2]
  }
  if (!chr %in% levels(chromosomes)) 
    stop("Chromosome not found. The input format should be chr#.")
  
  # Set vars
  dummy = subset(hm450, BiocGenerics::as.vector(chromosomes == chr))
  Windows = seq(from = from, to = to, by = winSize)
  Windows = cbind(Windows, append((Windows - 1)[-1], chrLengths[chr]))
  GR = GRanges(seqnames = chr, seqlengths = chrLengths[chr], ranges = IRanges(start = Windows[, 
                                                                                              1], end = Windows[, 2]))
  Ans = countOverlaps(GR, dummy)
  names(Ans) = paste(seqnames(GR), start(GR), end(GR), sep = ":")
  cat("Count complete for", chr, "\n")
  return(Ans)
}

.makeSeqCol = function(input, Gaps, pal, minVal = NULL, maxVal = NULL, ...) {
  minVal = ifelse(is.null(minVal), min(input), minVal)
  maxVal = ifelse(is.null(maxVal), max(input), maxVal)
  cuts = seq(minVal, maxVal, by = 10)
  groups = cut2(input, cuts = cuts, oneval = FALSE)
  outLegend = colorRampPalette(pal)(length(levels(groups)))
  names(outLegend) = levels(groups)
  outCol = outLegend[groups]
  Ans = list()
  Ans$colors = outCol
  Ans$legend = outLegend
  return(Ans)
}

.draw.acen = function(start, end, cytoband, ...) {
  rect(start, 0.25, end, 0.75, col = "maroon", ...)
  # segments(mean(c(start,end)),1,mean(c(start,end)),1.2)
  # text(mean(c(start,end)),1.1,mean(c(start,end)),1.1,labels=cytoband,srt=90,pos=3,cex=0.5)
}

.draw.giemsa = function(start, end, cytoband, ...) {
  rect(start, 0, end, 1, ...)
  segments(mean(c(start, end)), 1, mean(c(start, end)), 1.2)
  text(mean(c(start, end)), 1.1, mean(c(start, end)), 1.1, labels = cytoband, 
       srt = 90, pos = 3, cex = 0.5)
}

.drawIdeogram = function(x) {
  start = as.numeric(x[2])
  end = as.numeric(x[3])
  cytoband = x[4]
  giemsa = x[5]
  switch(giemsa, acen = .draw.acen(start, end, cytoband), 
         gneg = .draw.giemsa(start, end, cytoband, col = "white"), 
         gpos25 = .draw.giemsa(start, end, cytoband, col = "grey75"), 
         gpos50 = .draw.giemsa(start, end, cytoband, col = "grey60"), 
         gpos75 = .draw.giemsa(start, end, cytoband, col = "grey40"), 
         gpos100 = .draw.giemsa(start, end, cytoband, col = "grey15"), 
         gvar = .draw.giemsa(start, end, cytoband, col = "skyblue2"), 
         stalk = .draw.giemsa(start, end, cytoband, col = "tan"))
}