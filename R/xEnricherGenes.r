#' Function to conduct enrichment analysis given a list of genes and the ontology in query
#'
#' \code{xEnricherGenes} is supposed to conduct enrichment analysis given the input data and the ontology in query. It returns an object of class "eTerm". Enrichment analysis is based on either Fisher's exact test or Hypergeometric test. The test can respect the hierarchy of the ontology. Now it supports enrichment analysis using a wide variety of ontologies such as Gene Ontology and Phenotype Ontologies.
#'
#' @param data an input vector containing gene symbols
#' @param background a background vector containing gene symbols as the test background. If NULL, by default all annotatable are used as background
#' @param check.symbol.identity logical to indicate whether to match the input data/background via Synonyms for those unmatchable by official gene symbols. By default, it sets to false
#' @param ontology the ontology supported currently. By default, it is 'NA' to disable this option. Pre-built ontology and annotation data are detailed in \code{\link{xDefineOntology}}.
#' @param ontology.customised an object 'GS'. Higher priority over 'ontology' above. Required, otherwise it will return NULL
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the test statistic used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > fisher's exact test > binomial test. In other words, in terms of the calculated p-value, hypergeometric test < fisher's exact test < binomial test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to the annotatable. By default, it is NULL: if ontology.algorithm is not 'none', it is always TRUE; otherwise, it depends on the background (if not provided, it will be TRUE; otherwise FALSE). Surely, it can be explicitly stated
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note' below
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all genes in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is used to consider the children of a term (and subsequently all genes in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{term_info}: a matrix of nTerm X 4 containing snp/gene set information, where nTerm is the number of terms, and the 4 columns are "id" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{annotation}: a list of terms containing annotations, each term storing its annotations. Always, terms are identified by "id"}
#'  \item{\code{g}: an igraph object to represent DAG}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{background}: a vector containing the background data. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{overlap}: a list of overlapped snp/gene sets, each storing snps overlapped between a snp/gene set and the given input data (i.e. the snps of interest). Always, gene sets are identified by "id"}
#'  \item{\code{fc}: a vector containing fold changes}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{cross}: a matrix of nTerm X nTerm, with an on-diagnal cell for the overlapped-members observed in an individaul term, and off-diagnal cell for the overlapped-members shared betwene two terms}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The interpretation of the algorithms used to account for the hierarchy of the ontology is:
#' \itemize{
#' \item{"none": does not consider the ontology hierarchy at all.}
#' \item{"lea": computers the significance of a term in terms of the significance of its children at the maximum depth (e.g. 2). Precisely, once snps are already annotated to any children terms with a more signficance than itself, then all these snps are eliminated from the use for the recalculation of the signifance at that term. The final p-values takes the maximum of the original p-value and the recalculated p-value.}
#' \item{"elim": computers the significance of a term in terms of the significance of its all children. Precisely, once snps are already annotated to a signficantly enriched term under the cutoff of e.g. pvalue<1e-2, all these snps are eliminated from the ancestors of that term).}
#' \item{"pc": requires the significance of a term not only using the whole snps as background but also using snps annotated to all its direct parents/ancestors as background. The final p-value takes the maximum of both p-values in these two calculations.}
#' \item{"Notes": the order of the number of significant terms is: "none" > "lea" > "elim" > "pc".}
#' }
#' @export
#' @seealso \code{\link{xDefineOntology}}, \code{\link{xEnricher}}
#' @include xEnricherGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # Gene-based enrichment analysis using REACTOME pathways
#' # a) provide the input Genes of interest (eg 500 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 500))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # b) perform enrichment analysis
#' eTerm <- xEnricherGenes(data=data, ontology="MsigdbC2REACTOME", RData.location=RData.location)
#'
#' # c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' # d) save enrichment results to the file called 'REACTOME_enrichments.txt'
#' res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' output <- data.frame(term=rownames(res), res)
#' utils::write.table(output, file="REACTOME_enrichments.txt", sep="\t", row.names=FALSE)
#'
#' # e) barplot of significant enrichment results
#' gp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="adjp")
#' print(gp)
#'
#' # f) visualise the top 10 significant terms in the ontology hierarchy
#' # color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
#' xEnrichDAGplot(eTerm, top_num=10, displayBy="adjp", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25))
#' # color-code terms according to the z-scores
#' xEnrichDAGplot(eTerm, top_num=10, displayBy="zscore", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25))
#' 
#' # g) visualise the significant terms in the ontology hierarchy 
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- xRDataLoader(RData.customised='ig.REACTOME', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes=c("R-HSA-162582","R-HSA-168256"), mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#' xEnrichDAGplot(eTerm, top_num="auto", ig=ig, displayBy="adjp", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25))
#' }

xEnricherGenes <- function(data, background=NULL, check.symbol.identity=F, ontology=NA, ontology.customised=NULL, size.range=c(10,2000), min.overlap=5, which.distance=NULL, test=c("fisher","hypergeo","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    #ontology <- match.arg(ontology)
    ontology <- ontology[1]
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    
    ############
    if(length(data)==0){
    	return(NULL)
    }
    ############
    
    if (is.vector(data)){
        data <- unique(data)
    }else{
        warnings("The input data must be a vector.\n")
        return(NULL)
    }
    
    data <- as.character(data)
    
    #################################
 	aOnto <- xDefineOntology(ontology, ontology.customised=ontology.customised, verbose=verbose, RData.location=RData.location, guid=guid)
 	g <- aOnto$g
 	anno <- aOnto$anno
 	if(is.null(g)){
		warnings("There is no input for the ontology.\n")
        return(NULL)
	}
    #################################

    if(is.null(ontology.customised)){
		## convert gene symbol to entrz gene for both input data of interest and the input background (if given)
		if(verbose){
			now <- Sys.time()
			message(sprintf("Do gene mapping from Symbols to EntrezIDs (%s) ...", as.character(now)), appendLF=T)
		}
		data <- xSymbol2GeneID(data, check.symbol.identity=check.symbol.identity, verbose=verbose, RData.location=RData.location, guid=guid)
		data <- data[!is.na(data)]
		if(length(background)>0){
			background <- xSymbol2GeneID(background, check.symbol.identity=check.symbol.identity, verbose=verbose, RData.location=RData.location, guid=guid)
			background <- background[!is.na(background)]
		}
    }
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xEnricher' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    eTerm <- xEnricher(data=data, annotation=anno, g=g, background=background, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
	
	# replace EntrezGenes with gene symbols	
	if(is.null(ontology.customised) & class(eTerm)=="eTerm"){
		## load Enterz Gene information
		EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, guid=guid, verbose=verbose)
		allGeneID <- EG$gene_info$GeneID
		allSymbol <- as.vector(EG$gene_info$Symbol)
	
		## overlap
		overlap <- eTerm$overlap
		overlap_symbols <- lapply(overlap,function(x){
			ind <- match(x, allGeneID)
			allSymbol[ind]
		})
		eTerm$overlap <- overlap_symbols
		## data
		eTerm$data <- allSymbol[match(eTerm$data,allGeneID)]
		## background
		eTerm$background <- allSymbol[match(eTerm$background,allGeneID)]		
		
		## annotation
		annotation <- eTerm$annotation
		annotation_symbols <- lapply(annotation,function(x){
			ind <- match(x, allGeneID)
			allSymbol[ind]
		})
		eTerm$annotation <- annotation_symbols
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xEnricher' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xEnricherGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(eTerm)
}
