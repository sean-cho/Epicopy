#' Run Epicopy.
#' 
#' \code{epicopy} Master function to run epicopy from raw Illumina Methylation 450K data.
#' 
#' @description Parses the CNA result from \code{LRRtoCNA} and writes a
#'    segment file and marker file for input into GISTIC 2.0.
#' 
#' @param target_dir Target directory with idat files and sample sheet.
#' @param output_dir Output directory. Defaults to current directory. If FALSE,
#'    no GISTIC 2.0 files will be outputted, but can be manually obtained
#'    using the function \code{export_gistic}.
#' @param project_name Will be used as prefixes for any output.
#' @param Normals Character string identifying column in sample sheet to scan for
#'    normals. Defaults to NULL, using all included Epicopy normals. See 
#'    details for information on the normals and additional arguments.
#' @param sampNames Character string identifying column in sample sheet to use for
#'    sample names. If NULL, uses the default \code{read.450k.exp} values, which are
#'    the ChipNo_Position.
#' @param QN Logical. Quantile normalize post funnorm? Defaults to \code{FALSE}.
#' @param mode.bw If Ref = 'mode', uses this bw to estimate the mode. 
#'    Defaults to 0.1.
#' @param mode.method If Ref = 'mode', uses this method to estimate 
#'    the mode. Defaults to 'naive'.
#' @param normal.cnv Logical. Defaults to NULL which evaluates the Normal input.
#'    If Normal is user input normal tissues, then defaults to TRUE. Otherwise FALSE.
#' @param mean.center Logical. Perform mean centering of final data.
#' @param filterProbes Logical. Filter probe by TCGA defined 450K probes 
#'    adjacent to a known SNP.
#' @param retainedProbes Character. IF filterProbes is TRUE, denotes probes to keep. 
#'    If NULL, uses probes retained by TCGA.
#' @param ncores Integer. Number of cores to use for segmentation. Defaults to 1.
#' @param seg.undo.splits Undo.splits argument for \code{segment} from
#'    \code{DNAcopy}. Default: \code{'sdundo'}.
#' @param seg.undo.SD undo.SD argument for \code{segment}. Default: 2.
#'  @param filterbycount Logical. Recommended. Should the output segment file be 
#'    filtered for having at least \code{min_probes} number of probes in the 
#'    segment.
#' @param min_probes Number of probes to filter against.
#' @param verbose Logical. Verbose?
#' @param ... Passes argument to \code{getLRR}.
#' 
#' @details Epicopy reads the sample sheet provided by the user in the \code{target_dir},
#' imports the experiment, and returns a CNA file while writing gistic outputs at
#' the \code{output_dir}.
#' 
#' If the user has their own normals, specify the column name in the samplesheet
#' that contains the identifier 'normal' (not case sensitive). Else, they can
#' specify one of 'all' (default, NULL), 'breast', 'lung', or 'thyroid' to use Epicopy 
#' included normals. The last alternative is NA which uses all the samples as 
#' references, under the assumption that the \code{mode}/\code{median} of all the
#' samples are at 2n. The last is recommended either in samples with low CNV, or
#' when there are many samples. See also \code{\link[Epicopy]{getLRR}}.
#' 
#' @import BiocGenerics Hmisc GenomeInfoDb GenomicRanges IRanges modeest minfi EpicopyData
#' @export 
#' 
#' @seealso \code{\link[Epicopy]{getLRR}} \code{\link[Epicopy]{LRRtoCNA}} 
#'    \code{\link[Epicopy]{export_gistic}} \code{\link[minfi]{RGChannelSet-class}}
#'    \code{\link[minfi]{read.450k.sheet}} \code{\link[minfi]{read.450k.exp}}

epicopy <- function(target_dir, output_dir = NULL, project_name = NULL,
                    Normals = NULL, sampNames = NULL, QN = FALSE,
                    Ref = 'mode', mode.bw = 0.1, mode.method = "naive", 
                    normal.cnv = NULL, mean.center = TRUE, 
                    filterProbes = FALSE, retainedProbes = NULL, 
                    ncores = 1L, seg.undo.splits = 'sdundo', seg.undo.SD = 2,
                    filterbycount = TRUE, min_probes = 50, verbose = TRUE, ...){
  
  if(is.null(output_dir)){
    cat('Using current directory as output directory.\n')
    output_dir <- '.'
  }
  
  target_sheet <- read.450k.sheet(target_dir)

  if(!is.null(Normals)){
    if(!Normals %in% c('breast', 'lung', 'thyroid', 'all')){
      if(!inherits(Normals, 'character')) stop('Normals input not recognized. Please review and see ?epicopy.\n')
      cat('Looking for normal annotation under', Normals, 'column in samplesheet.\n')
      normal_index <- tolower(target_sheet[,Normals])
      if(!any(normal_index %in% 'normal'))
        stop('No normal sample found in annotation.\n')
      nnormals <- table(normal_index %in% 'normal')['TRUE']
      cat('Found', nnormals, 'normal samples.\n')
      Normals <- normal_index %in% 'normal'
      if(is.null(normal.cnv)) normal.cnv <- TRUE
    } else {
      if(is.null(normal.cnv)) normal.cnv <- FALSE
    }
  } else {
    if(is.null(normal.cnv)) normal.cnv <- FALSE
  }
  
  rgset <- read.450k.exp(targets = target_sheet, 
                         verbose = verbose)
  
  lrr <- getLRR(rgSet = rgset, Normals = Normals, sampNames = sampNames, QN = QN,
         Ref = Ref, mode.bw = mode.bw, mode.method = mode.method, 
         normal.cnv = normal.cnv, mean.center = mean.center, 
         filterProbes = filterProbes, ...)
  
  cna <- LRRtoCNA(lrr, ncores = ncores, seg.undo.splits = seg.undo.splits,
                  seg.undo.SD = seg.undo.SD)

  if(!is.null(output_dir)){
    if(output_dir == FALSE){
      return(cna)
    }
  }

  export_gistic(input = cna, output_dir = output_dir, filterbycount = filterbycount,
                min_probes = min_probes, segfile_name = project_name,
                markerfile_name = project_name)
  
  return(cna)
}