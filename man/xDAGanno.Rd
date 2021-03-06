% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xDAGanno.r
\name{xDAGanno}
\alias{xDAGanno}
\title{Function to generate a subgraph of a direct acyclic graph (DAG) induced by the input annotation data}
\usage{
xDAGanno(
g,
annotation,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = TRUE,
verbose = TRUE
)
}
\arguments{
\item{g}{an object of class "igraph" to represent DAG}

\item{annotation}{the vertices/nodes for which annotation data are
provided. It can be a sparse Matrix of class "dgCMatrix" (with
variants/genes as rows and terms as columns), or a list of nodes/terms
each containing annotation data, or an object of class 'GS' (basically
a list for each node/term with annotation data)}

\item{path.mode}{the mode of paths induced by vertices/nodes with input
annotation data. It can be "all_paths" for all possible paths to the
root, "shortest_paths" for only one path to the root (for each node in
query), "all_shortest_paths" for all shortest paths to the root (i.e.
for each node, find all shortest paths with the equal lengths)}

\item{true.path.rule}{logical to indicate whether the true-path rule
should be applied to propagate annotations. By default, it sets to
true}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to true for display}
}
\value{
\itemize{
\item{\code{subg}: an induced subgraph, an object of class "igraph". In
addition to the original attributes to nodes and edges, the return
subgraph is also appended by two node attributes: 1) "anno" containing
a list of variants/genes either as original annotations (and inherited
annotations; 2) "IC" standing for information content defined as
negative 10-based log-transformed frequency of variants/genes annotated
to that term.}
}
}
\description{
\code{xDAGanno} is supposed to produce a subgraph induced by the input
annotation data, given a direct acyclic graph (DAG; an ontology). The
input is a graph of "igraph", a list of the vertices containing
annotation data, and the mode defining the paths to the root of DAG.
The induced subgraph contains vertices (with annotation data) and their
ancestors along with the defined paths to the root of DAG. The
annotations at these vertices (including their ancestors) can also be
updated according to the true-path rule: those annotated to a term
should also be annotated by its all ancestor terms.
}
\note{
For the mode "shortest_paths", the induced subgraph is the most
concise, and thus informative for visualisation when there are many
nodes in query, while the mode "all_paths" results in the complete
subgraph.
}
\examples{
\dontrun{
# 1) SNP-based ontology
# 1a) ig.EF (an object of class "igraph" storing as a directed graph)
g <- xRDataLoader('ig.EF')

# 1b) load GWAS SNPs annotated by EF (an object of class "dgCMatrix" storing a spare matrix)
anno <- xRDataLoader(RData='GWAS2EF')

# 1c) prepare for annotation data
# randomly select 5 terms/vertices (and their annotation data)
annotation <- anno[, sample(1:dim(anno)[2],5)]

# 1d) obtain the induced subgraph according to the input annotation data
# based on shortest paths (i.e. the most concise subgraph induced)
dag <- xDAGanno(g, annotation, path.mode="shortest_paths",
verbose=TRUE)

# 1e) color-code nodes/terms according to the number of annotations
data <- sapply(V(dag)$anno, length)
names(data) <- V(dag)$name
dnet::visDAG(g=dag, data=data, node.info="both")

####################
# Below is for those SNPs annotated by the term called 'ankylosing spondylitis'
# The steps 1a) and 1b) are the same as above
# 1c') prepare for annotation data
# select a term 'ankylosing spondylitis'
terms <- V(g)$term_id[grep('ankylosing spondylitis',V(g)$term_name,
perl=TRUE)]
ind <- which(colnames(anno) \%in\% terms)
annotation <- lapply(ind, function(x){names(which(anno[,x]!=0))})
names(annotation) <- colnames(anno)[ind]

# 1d') obtain the induced subgraph according to the input annotation data
# based on all possible paths (i.e. the complete subgraph induced)
dag <- xDAGanno(g, annotation, path.mode="all_paths", verbose=TRUE)

# 1e') color-code nodes/terms according to the number of annotations
data <- sapply(V(dag)$anno, length)
names(data) <- V(dag)$name
dnet::visDAG(g=dag, data=data, node.info="both")

###########################################################
# 2) Gene-based ontology
# 2a) ig.MP (an object of class "igraph" storing as a directed graph)
g <- xRDataLoader('ig.MP')

# 2b) load human genes annotated by MP (an object of class "GS" containing the 'gs' component)
GS <- xRDataLoader(RData='org.Hs.egMP')
anno <- GS$gs # notes: This is a list

# 2c) prepare for annotation data
# randomly select 5 terms/vertices (and their annotation data)
annotation <- anno[sample(1:length(anno),5)]

# 2d) obtain the induced subgraph according to the input annotation data
# based on shortest paths (i.e. the most concise subgraph induced) 
# but without applying true-path rule
dag <- xDAGanno(g, annotation, path.mode="shortest_paths",
true.path.rule=TRUE, verbose=TRUE)

# 2e) color-code nodes/terms according to the number of annotations
data <- sapply(V(dag)$anno, length)
names(data) <- V(dag)$name
dnet::visDAG(g=dag, data=data, node.info="both")
}
}
\seealso{
\code{\link{xRDataLoader}}
}
