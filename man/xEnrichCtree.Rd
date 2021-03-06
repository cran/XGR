% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xEnrichCtree.r
\name{xEnrichCtree}
\alias{xEnrichCtree}
\title{Function to visualise enrichment results using a tree-like circular plot}
\usage{
xEnrichCtree(
eTerm,
ig,
FDR.cutoff = NULL,
node.color = c("zscore", "adjp", "or", "nOverlap"),
colormap = "brewer.Reds",
zlim = NULL,
node.size = c("adjp", "zscore", "or", "nOverlap"),
slim = NULL,
node.size.range = c(0.5, 4.5),
group.gap = 0.08,
group.color = "lightblue",
group.size = 0.2,
group.label.size = 2,
group.label.color = "black",
legend.direction = c("auto", "horizontal", "vertical"),
leave.label.orientation = c("inwards", "outwards"),
...
)
}
\arguments{
\item{eTerm}{an object of class "eTerm" or "ls_eTerm". Alterntively, it
can be a data frame having all these columns
('name','adjp','or','zscore','nOverlap'; 'group' optionally)}

\item{ig}{an object of class "igraph" with node attribute 'name'. It
could be a 'phylo' object converted to. Note: the leave labels would be
the node attribute 'name' unless the node attribute 'label' is
explicitely provided}

\item{FDR.cutoff}{FDR cutoff used to show the significant terms only.
By default, it is set to NULL; useful when nodes sized by FDR}

\item{node.color}{which statistics will be used for node coloring. It
can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and
"zscore" for enrichment z-score}

\item{colormap}{short name for the colormap. It can be one of "jet"
(jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red
colormap), "wyr" (white-yellow-red colormap), "br" (black-red
colormap), "yr" (yellow-red colormap), "wb" (white-black colormap),
"rainbow" (rainbow colormap, that is,
red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2
default color palette). Alternatively, any hyphen-separated HTML color
names, e.g. "lightyellow-orange" (by default), "blue-black-yellow",
"royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of
standard color names can be found in
\url{http://html-color-codes.info/color-names}}

\item{zlim}{the minimum and maximum values for which colors should be
plotted}

\item{node.size}{which statistics will be used for node size. It can be
"or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore"
for enrichment z-score}

\item{slim}{the minimum and maximum values for which sizes should be
plotted}

\item{node.size.range}{the range of actual node size}

\item{group.gap}{the gap between group circles. Only works when
multiple groups provided}

\item{group.color}{the color of group circles. Only works when multiple
groups provided}

\item{group.size}{the line width of group circles. Only works when
multiple groups provided}

\item{group.label.size}{the size of group circle labelling. Always a
sequential integer located at the top middle. Only works when multiple
groups provided}

\item{group.label.color}{the color of group circle labelling. Only
works when multiple groups provided}

\item{legend.direction}{the legend guide direction. It can be
"horizontal" (useful for many groups with lengthy labelling),
"vertical" and "auto" ("vertical" when multiple groups provided;
otherwise "horizontal")}

\item{leave.label.orientation}{the leave label orientation. It can be
"outwards" and "inwards"}

\item{...}{additional graphic parameters used in xCtree}
}
\value{
a ggplot2 object appended with 'ig', 'data' which should contain
columns 'x','y', 'leaf' (T/F), 'name' (the same as V(ig)$name), 'tipid'
(tip id), 'label' (if not given in ig, a 'name' varient), 'angle' and
'hjust' (assist in leave label orientation), and 'data_enrichment'
(enrichment results for tips)
}
\description{
\code{xEnrichCtree} is supposed to visualise enrichment results using a
tree-like circular plot.
}
\note{
none
}
\examples{
\dontrun{
# Load the XGR package and specify the location of built-in data
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

# load the atlas of AA pathways
AA.template <- xRDataLoader("AA.template",
RData.location=RData.location)
# consensus tree
ig <- AA.template$consensus$ig

# enrichment analysis using AA pathways
input <- xRDataLoader('Haploid_regulators_all',
RData.location=RData.location)
data <- subset(input, Phenotype=="AKT")
genes <- data$Gene[data$FDR<0.05]
background <- data$Gene
eTerm <- xEnricherGenes(genes, background=background, ontology="AA",
min.overlap=5, test="fisher", RData.location=RData.location)

# circular visualisation of enriched AA pathways
gp <- xEnrichCtree(eTerm, ig)

###############################
# advanced use: multiple groups
# enrichment analysis using AA pathways
Haploid <- subset(input, FDR<0.05)
ls_group <- split(x=Haploid$Gene, f=Haploid$Phenotype)
background <- unique(input$Gene)
ls_eTerm <- xEnricherGenesAdv(ls_group, background=background,
ontologies="AA", test="fisher", min.overlap=5,
RData.location=RData.location)

# circular visualisation of enriched AA pathways
gp <- xEnrichCtree(ls_eTerm, ig)
}
}
\seealso{
\code{\link{xCtree}}
}
