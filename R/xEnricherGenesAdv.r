#' Function to conduct enrichment analysis given a list of gene sets and a list of ontologies
#'
#' \code{xEnricherGenesAdv} is supposed to conduct enrichment analysis given a list of gene sets and a list of ontologies. It is an advanced version of \code{xEnricherGenes}, returning an object of the class 'ls_eTerm'.
#'
#' @param list_vec an input vector containing gene symbols. Alternatively it can be a list of vectors, representing multiple groups of genes
#' @param background a background vector containing gene symbols as the test background. If NULL, by default all annotatable are used as background
#' @param check.symbol.identity logical to indicate whether to match the input data/background via Synonyms for those unmatchable by official gene symbols. By default, it sets to false
#' @param ontologies the ontologies supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PSG" for phylostratigraphy (phylostratific age), "PS" for sTOL-based phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for SCOP domain superfamilies, "Pfam" for Pfam domain families, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, "EF" for Experimental Factor Ontology (used to annotate GWAS Catalog genes), Drug-Gene Interaction database ("DGIdb") for druggable categories, tissue-specific eQTL-containing genes from GTEx ("GTExV4", "GTExV6p" and "GTExV7"), crowd extracted expression of differential signatures from CREEDS ("CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP" and "CreedsGeneDN"), KEGG pathways (including 'KEGG' for all, 'KEGGmetabolism' for 'Metabolism' pathways, 'KEGGgenetic' for 'Genetic Information Processing' pathways, 'KEGGenvironmental' for 'Environmental Information Processing' pathways, 'KEGGcellular' for 'Cellular Processes' pathways, 'KEGGorganismal' for 'Organismal Systems' pathways, and 'KEGGdisease' for 'Human Diseases' pathways), 'REACTOME' for REACTOME pathways or 'REACTOME_x' for its sub-ontologies (where x can be 'CellCellCommunication', 'CellCycle', 'CellularResponsesToExternalStimuli', 'ChromatinOrganization', 'CircadianClock', 'DevelopmentalBiology', 'DigestionAndAbsorption', 'Disease', 'DNARepair', 'DNAReplication', 'ExtracellularMatrixOrganization', 'GeneExpression(Transcription)', 'Hemostasis', 'ImmuneSystem', 'Metabolism', 'MetabolismOfProteins', 'MetabolismOfRNA', 'Mitophagy', 'MuscleContraction', 'NeuronalSystem', 'OrganelleBiogenesisAndMaintenance', 'ProgrammedCellDeath', 'Reproduction', 'SignalTransduction', 'TransportOfSmallMolecules', 'VesicleMediatedTransport'), and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
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
#' @param plot logical to indicate whether heatmap plot is drawn
#' @param fdr.cutoff fdr cutoff used to declare the significant terms. By default, it is set to 0.05. This option only works when setting plot (see above) is TRUE
#' @param displayBy which statistics will be used for drawing heatmap. It can be "fc" for enrichment fold change, "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score (by default), "or" for odds ratio. This option only works when setting plot (see above) is TRUE
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "ls_eTerm", a list with following components:
#' \itemize{
#'  \item{\code{df}: a data frame of n x 12, where the 12 columns are "group" (the input group names), "ontology" (input ontologies), "id" (term ID), "name" (term name), "nAnno" (number in members annotated by a term), "nOverlap" (number in overlaps), "fc" (enrichment fold changes), "zscore" (enrichment z-score), "pvalue" (nominal p value), "adjp" (adjusted p value (FDR)), "or" (odds ratio), "CIl" (lower bound confidence interval for the odds ratio), "CIu" (upper bound confidence interval for the odds ratio), "distance" (term distance or other information), "members" (members (represented as Gene Symbols) in overlaps)}
#'  \item{\code{mat}: NULL if the plot is not drawn; otherwise, a matrix of term names X groups with numeric values for the signficant enrichment, NA for the insignificant ones}
#'  \item{\code{gp}: NULL if the plot is not drawn; otherwise, a 'ggplot' object}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xEnricherGenes}}, \code{\link{xEnrichViewer}}, \code{\link{xHeatmap}}
#' @include xEnricherGenesAdv.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # Gene-based enrichment analysis using ontologies (REACTOME and GOMF)
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # b) perform enrichment analysis
#' ls_eTerm <- xEnricherGenesAdv(data, ontologies=c("REACTOME","GOMF"), RData.location=RData.location)
#' ls_eTerm
#' ## forest plot of enrichment results
#' gp <- xEnrichForest(ls_eTerm, top_num=10)
#' }

xEnricherGenesAdv <- function(list_vec, background=NULL, check.symbol.identity=F, ontologies=c("GOBP","GOMF","GOCC","PSG","PS","PS2","SF","Pfam","DO","HPPA","HPMI","HPCM","HPMA","MP", "EF", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6p", "GTExV7", "CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN", "KEGG","KEGGmetabolism","KEGGgenetic","KEGGenvironmental","KEGGcellular","KEGGorganismal","KEGGdisease", "REACTOME", "REACTOME_ImmuneSystem", "REACTOME_SignalTransduction", "CGL"), size.range=c(10,2000), min.overlap=3, which.distance=NULL, test=c("hypergeo","fisher","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, verbose=T, silent=FALSE, plot=TRUE, fdr.cutoff=0.05, displayBy=c("zscore","fdr","pvalue","fc","or"), RData.location="http://galahad.well.ox.ac.uk/bigdata")
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
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    displayBy <- match.arg(displayBy)
    
    ############
    if(length(list_vec)==0){
    	return(NULL)
    }
    ############
    if(is.vector(list_vec) & class(list_vec)!="list"){
    	list_vec <- list(list_vec)
	}else if(class(list_vec)=="list"){
		## Remove null elements in a list
		list_vec <- base::Filter(base::Negate(is.null), list_vec)
		if(length(list_vec)==0){
			return(NULL)
		}	
    }else{
        stop("The input data must be a vector or a list of vectors.\n")
    }
    
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('G', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    ls_df <- lapply(1:length(list_vec), function(i){
		
		if(verbose){
			message(sprintf("Analysing group %d ('%s') (%s) ...", i, names(list_vec)[i], as.character(Sys.time())), appendLF=T)
		}
		data <- list_vec[[i]]
    	
    	ls_df <- lapply(1:length(ontologies), function(j){
			if(verbose){
				message(sprintf("\tontology %d ('%s') (%s) ...", j, ontologies[j], as.character(Sys.time())), appendLF=T)
			}
			ontology <- ontologies[j]
			
    		eTerm <- xEnricherGenes(data=data, background=background, check.symbol.identity=check.symbol.identity, ontology=ontology, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose, silent=TRUE, RData.location=RData.location)
			df <- xEnrichViewer(eTerm, top_num="all", sortBy="or", details=TRUE)
			
			if(is.null(df)){
				return(NULL)
			}else{
				cbind(group=rep(names(list_vec)[i],nrow(df)), ontology=rep(ontology,nrow(df)), id=rownames(df), df, stringsAsFactors=F)
			}
		})
		df <- do.call(rbind, ls_df)
	})
    df_all <- do.call(rbind, ls_df)
    
    ## heatmap view
    if(plot & !is.null(df_all)){

    	adjp <- NULL
    	
    	gp <- NULL
    	mat <- NULL
    	
    	ls_df <- split(x=df_all[,-12], f=df_all$ontology)
    	#######
    	## keep the same order for ontologies as input
    	ls_df <- ls_df[unique(df_all$ontology)]
    	#######
		ls_mat <- lapply(1:length(ls_df), function(i){
		
			df <- ls_df[[i]]
			ind <- which(df$adjp < fdr.cutoff)
			if(length(ind)>=1){
				df <- as.data.frame(df %>% dplyr::filter(adjp < fdr.cutoff))
				
				if(displayBy=='fdr'){
					mat <- as.matrix(xSparseMatrix(df[,c('name','group','adjp')], rows=unique(df$name), columns=names(list_vec)))
					mat[mat==0] <- NA
					mat <- -log10(mat)
				}else if(displayBy=='pvalue'){
					mat <- as.matrix(xSparseMatrix(df[,c('name','group','pvalue')], rows=unique(df$name), columns=names(list_vec)))
					mat[mat==0] <- NA
					mat <- -log10(mat)
				}else if(displayBy=='zscore'){
					mat <- as.matrix(xSparseMatrix(df[,c('name','group','zscore')], rows=unique(df$name), columns=names(list_vec)))
					mat[mat==0] <- NA
				}else if(displayBy=='fc'){
					mat <- as.matrix(xSparseMatrix(df[,c('name','group','fc')], rows=unique(df$name), columns=names(list_vec)))
					mat[mat==0] <- NA
					mat <- log2(mat)
				}else if(displayBy=='or'){
					mat <- as.matrix(xSparseMatrix(df[,c('name','group','or')], rows=unique(df$name), columns=names(list_vec)))
					mat[mat==0] <- NA
					mat <- log2(mat)
				}
				
				if(nrow(mat)==1){
					df_mat <- mat
				}else{
					
					## order by the length of names
					rname_ordered <- rownames(mat)[order(-nchar(rownames(mat)))]
					## order by the evolutionary ages
					if(names(ls_df)[i]=='PS2' || names(ls_df)[i]=='PSG'){
						df_tmp <- unique(df[,c('id','name')])
						df_tmp <- df_tmp[with(df_tmp, order(as.numeric(df_tmp$id))),]
						rname_ordered <- df_tmp$name
					}
					
					ind <- match(rname_ordered, rownames(mat))
					df_mat <- as.matrix(mat[ind,], ncol=ncol(mat))
					colnames(df_mat) <- colnames(mat)
					
					colnames(df_mat) <- colnames(mat)
				}
				return(df_mat)
					
			}else{
				return(NULL)
			}
			
		})
		mat <- do.call(rbind, ls_mat)
		
		if(!is.null(mat)){
			if(displayBy=='fdr' | displayBy=='pvalue'){
				colormap <- 'grey100-darkorange'
				zlim <- c(0, ceiling(max(mat[!is.na(mat)])))
				
				if(displayBy=='fdr'){
					legend.title <- expression(-log[10]("FDR"))
				}else if(displayBy=='pvalue'){
					legend.title <- expression(-log[10]("p-value"))
				}
				
			}else if(displayBy=='fc' | displayBy=='zscore' | displayBy=='or'){
				tmp_max <- ceiling(max(mat[!is.na(mat)]))
				tmp_min <- floor(min(mat[!is.na(mat)]))
				if(tmp_max>0 & tmp_min<0){
					colormap <- 'deepskyblue-grey100-darkorange'	
					tmp <- max(tmp_max, abs(tmp_min))
					zlim <- c(-tmp, tmp)
				}else if(tmp_max<=0){
					colormap <- 'deepskyblue-grey100'	
					zlim <- c(tmp_min, 0)
				}else if(tmp_min>=0){
					colormap <- 'grey100-darkorange'	
					zlim <- c(0, tmp_max)
				}
				
				if(displayBy=='fc'){
					legend.title <- expression(log[2]("FC"))
				}else if(displayBy=='zscore'){	
					legend.title <- ("Z-score")
				}else if(displayBy=='or'){	
					legend.title <- expression(log[2]("OR"))
				}
			}
		
			gp <- xHeatmap(mat, reorder="none", colormap=colormap, ncolors=64, zlim=zlim, legend.title=legend.title, barwidth=0.4, x.rotate=60, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent',barheight=max(3,min(5,nrow(mat))))
			gp <- gp + theme(legend.title=element_text(size=8))
		}
		
    }else{
    	mat <- NULL
    	gp <- NULL
    }
    
    ls_eTerm <- list(df = df_all,
    			   mat = mat,
    			   gp = gp
                 )
    class(ls_eTerm) <- "ls_eTerm"
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(ls_eTerm)
}
