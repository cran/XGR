% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xSNP2eGenes.r
\name{xSNP2eGenes}
\alias{xSNP2eGenes}
\title{Function to define eQTL genes given a list of SNPs or a customised eQTL mapping data}
\usage{
xSNP2eGenes(
data,
include.eQTL = NA,
eQTL.customised = NULL,
cdf.function = c("empirical", "exponential"),
plot = FALSE,
verbose = TRUE,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL
)
}
\arguments{
\item{data}{an input vector containing SNPs. SNPs should be provided as
dbSNP ID (ie starting with rs). Alternatively, they can be in the
format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for
example, 'chr16:28525386'}

\item{include.eQTL}{the eQTL supported currently. By default, it is
'NA' to disable this option. Pre-built eQTL datasets are detailed in
\code{\link{xDefineEQTL}}}

\item{eQTL.customised}{a user-input matrix or data frame with 4
columns: 1st column for SNPs/eQTLs, 2nd column for Genes, 3rd for eQTL
mapping significance level (p-values or FDR), and 4th for contexts
(required even though only one context is input). Alternatively, it can
be a file containing these 4 columns. It is designed to allow the user
analysing their eQTL data. This customisation (if provided) will
populate built-in eQTL data}

\item{cdf.function}{a character specifying a Cumulative Distribution
Function (cdf). It can be one of 'exponential' based on exponential
cdf, 'empirical' for empirical cdf}

\item{plot}{logical to indicate whether the histogram plot (plus
density or CDF plot) should be drawn. By default, it sets to false for
no plotting}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to true for display}

\item{RData.location}{the characters to tell the location of built-in
RData files. See \code{\link{xRDataLoader}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF
project. See \code{\link{xRDataLoader}} for details}
}
\value{
a data frame with following columns:
\itemize{
\item{\code{Gene}: eQTL-containing genes}
\item{\code{SNP}: eQTLs}
\item{\code{Sig}: the eQTL mapping significant level (the
best/minimum)}
\item{\code{Weight}: the eQTL weight}
}
}
\description{
\code{xSNP2eGenes} is supposed to define eQTL genes given a list of
SNPs or a customised eQTL mapping data. The eQTL weight is calcualted
as Cumulative Distribution Function of negative log-transformed
eQTL-reported signficance level.
}
\note{
none
}
\examples{
\dontrun{
# Load the library
library(XGR)
}

RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
\dontrun{
# a) provide the SNPs with the significance info
data(ImmunoBase)
gr <- ImmunoBase$AS$variants
AS <- as.data.frame(GenomicRanges::mcols(gr)[, c('Variant','Pvalue')])

# b) define eQTL genes
df_eGenes <- xSNP2eGenes(data=AS[,1], include.eQTL="JKscience_TS2A",
RData.location=RData.location)
}
}
\seealso{
\code{\link{xRDataLoader}}
}
