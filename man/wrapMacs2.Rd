% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapmacs2.R
\name{wrapMacs2}
\alias{wrapMacs2}
\title{wrapMacs2}
\usage{
wrapMacs2(
  sce,
  outdir,
  genome = "mm",
  extsize = 200,
  shift = -(extsize/2),
  broadCall = "--broad",
  bampath = NULL
)
}
\arguments{
\item{sce}{A SingleCellExperiment object to populate with macs2 peaks}

\item{outdir}{character for the path where to store the peak files}

\item{genome}{the genome to use (default is `mm`)}

\item{extsize}{number of bp to extend the reads (see macs2 callpeak for more
info)}

\item{shift}{it is used by macs2 to cut the ends (5') towards 5'->3'
direction then apply EXTSIZE to extend them to fragments.}

\item{broadCall}{uses the --broad option, for broad peaks.}

\item{bampath}{the path to the bamfiles, in case of NULL it attempts to
retrieve the path from the `metadata(sce)$ct_bams` (default is NULL)}
}
\value{
a SingleCellExperiment populated with the new peaks
}
\description{
wrapMacs2
}
\examples{
TBD
}
