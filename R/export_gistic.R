#' Plot coverage of 450K across hg19
#' 
#' \code{export_gistic} Exports gistic inputs.
#' 
#' @description Parses the CNA result from \code{LRRtoCNA} and writes a
#'    segment file and marker file for input into GISTIC 2.0.
#' 
#' @param input Result from \code{LRRtoCNA}.
#' @param output Output directory.
#' @param filterbycount Logical. Recommended. Should the output segment file be 
#'    filtered for having at least \code{min_probes} number of probes in the 
#'    segment.
#' @param min_probes Number of probes to filter against.
#' @param segfile_name Name of the output segment file. NULL prints segmentation_output.
#' @param marker_file Name of the output marker file. NULL prints marker_file.
#' 
#' @details Argument \code{min_probes} restricts segments to a certain number of
#'    probes or more if \code{filterbycount} is \code{TRUE}. Writes two output files,
#'    segmentation_output and marker_file which are the .seg and marker file
#'    inputs required for GISTIC 2.0.
#'    
#' @examples
#' 
#' # Load included sample segmentation data
#' data(data_vignette)
#' # Export gistic files (defaults to local directory)
#' export_gistic(epi_seg, filterbycount = TRUE, min_probes = 50)
#'
#' @import BiocGenerics Hmisc GenomeInfoDb GenomicRanges IRanges
#' @export export_gistic
#' 
#' @seealso \code{\link[minfi]{preprocessFunnorm}} 
#'    \code{\link[minfi]{RGChannelSet-class}}

export_gistic <- function(input, output_dir, filterbycount = TRUE,
                          min_probes = 50, segfile_name = NULL,
                          markerfile_name = NULL){
  probeset <- Epicopy:::hm450[rownames(input$data),]
  marker_file <- data.frame(Probe_ID = names(probeset),
                            Chromosome = as.character(seqnames(probeset)),
                            Position = probeset$probeTarget)
  
  segfile <- input$output
  colnames(segfile) <- c('ID', 'Chromosome', 'Start', 'End', 'Probe_count', 'Segment_Mean')
  segfile$Chromosome <- paste0('chr', segfile$Chromosome)
  
  if(filterbycount) segfile <- subset(segfile, Probe_count >= min_probes)
  
  # Write .seg file
  if(is.null(segfile_name)) segfile_name <- 'segmentation_output.txt'
  segfile_name <- paste0(output_dir, '/', segfile_name)
  write.table(segfile, file = segfile_name,
              quote = FALSE, sep = '\t', row.names = FALSE)
  cat('Segmentation file written.\n')
  
  # Write marker file
  if(is.null(markerfile_name)) markerfile_name <- 'epicopy_markers.txt'
  markerfile_name <- paste0(output_dir, '/', markerfile_name)
  if(markerfile_name == segfile_name) {
    cat('Marker file and seg file names are equal. 
        Changing marker file to default name.\n')
    markerfile_name <- 'epicopy_markers.txt'
    markerfile_name <- paste0(output_dir, markerfile_name)
  }
  write.table(marker_file, file = markerfile_name,
              quote = FALSE, sep = '\t', row.names = FALSE)
  cat('Marker file written.\n')
}

