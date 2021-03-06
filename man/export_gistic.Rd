% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_gistic.R
\name{export_gistic}
\alias{export_gistic}
\title{Export GISTIC-compatible results}
\usage{
export_gistic(input, output_dir = NULL, filterbycount = TRUE,
  min_probes = 50, segfile_name = NULL, markerfile_name = NULL)
}
\arguments{
\item{input}{Result from \code{LRRtoCNA}.}

\item{output_dir}{Output directory.}

\item{filterbycount}{Logical. Recommended. Should the output segment file be 
filtered for having at least \code{min_probes} number of probes in the 
segment.}

\item{min_probes}{Number of probes to filter against.}

\item{segfile_name}{Name of the output segment file. NULL prints segmentation_output.}

\item{markerfile_name}{Name of the output marker file. NULL prints marker_file.}
}
\value{
Exports segmented data into GISTIC2 formatted input files.
}
\description{
Parses the CNA result from \code{LRRtoCNA} and writes a
   segment file and marker file for input into GISTIC 2.0.
}
\details{
\code{export_gistic} Exports gistic inputs.

Argument \code{min_probes} restricts segments to a certain number of
   probes or more if \code{filterbycount} is \code{TRUE}. Writes two output files,
   segmentation_output and marker_file which are the .seg and marker file
   inputs required for GISTIC 2.0.
}
\examples{

# Load included sample segmentation data
data(data_vignette)
# Export gistic files (defaults to local directory)
export_gistic(epi_seg, filterbycount = TRUE, min_probes = 50)

}
\seealso{
\code{\link[minfi]{preprocessFunnorm}} 
   \code{\link[minfi]{RGChannelSet-class}}
}
