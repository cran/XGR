% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xLiftOver.r
\name{xLiftOver}
\alias{xLiftOver}
\title{Function to lift genomic intervals from one genome build to another.}
\usage{
xLiftOver(
data.file,
format.file = c("data.frame", "bed", "chr:start-end", "GRanges"),
build.conversion = c(NA, "hg38.to.hg19", "hg19.to.hg38",
"hg19.to.hg18",
"hg18.to.hg38", "hg18.to.hg19"),
merged = T,
verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL
)
}
\arguments{
\item{data.file}{an input data file, containing a list of genomic
regions to test. If the input file is formatted as a 'data.frame'
(specified by the parameter 'format.file' below), the first three
columns correspond to the chromosome (1st column), the starting
chromosome position (2nd column), and the ending chromosome position
(3rd column). If the format is indicated as 'bed' (browser extensible
data), the same as 'data.frame' format but the position is 0-based
offset from chromomose position. If the genomic regions provided are
not ranged but only the single position, the ending chromosome position
(3rd column) is allowed not to be provided. If the format is indicated
as "chr:start-end", instead of using the first 3 columns, only the
first column will be used and processed. If the file also contains
other columns, these additional columns will be ignored. Alternatively,
the input file can be the content itself assuming that input file has
been read. Note: the file should use the tab delimiter as the field
separator between columns}

\item{format.file}{the format for input files. It can be one of
"data.frame", "chr:start-end", "bed"}

\item{build.conversion}{the conversion from one genome build to
another. The conversions supported are "hg38.to.hg19", "hg19.to.hg38",
"hg19.to.hg18", "hg18.to.hg38" and "hg18.to.hg19". By default it is NA,
forcing the user to specify the corrent one.}

\item{merged}{logical to indicate whether multiple ranges should be
merged into the one per a range in query. By default, it sets to true}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to false for no display}

\item{RData.location}{the characters to tell the location of built-in
RData files. See \code{\link{xRDataLoader}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF
project. See \code{\link{xRDataLoader}} for details}
}
\value{
an GR oject storing converted genomic intervals.
}
\description{
\code{xLiftOver} is supposed to lift genomic intervals from one genome
build to another. Supported are the conversions between genome builds
'hg38' (GRCh38), 'hg19' (GRCh37) and 'h18'.
}
\examples{
\dontrun{
# Load the XGR package and specify the location of built-in data
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"

# Provide UCSC genes (hg19)
UCSC_genes <- xRDataLoader(RData.customised='UCSC_genes',
RData.location=RData.location)
UCSC_genes

# Lift over to hg38
gr <- xLiftOver(UCSC_genes, format.file="GRanges",
build.conversion="hg19.to.hg38", RData.location=RData.location)
gr
}
}
\seealso{
\code{\link{xLiftOver}}
}
