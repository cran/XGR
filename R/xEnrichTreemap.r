#' Function to visualise enrichment results using a treemap
#'
#' \code{xEnrichTreemap} is supposed to visualise enrichment results using a treemap. The area is proportional to odds ratio, colored by the significance level. It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns (named as 'group','ontology','name','zscore','adjp','or','CIl','CIu','nOverlap','members_Overlap')
#' @param top_num the number of the top terms (sorted according to OR). If it is 'auto', only the significant terms (see below FDR.cutoff) will be displayed
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param CI.one logical to indicate whether to allow the inclusion of one in CI. By default, it is TURE (allowed)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the -log10(FDR)
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param wrap.width a positive integer specifying wrap width of name
#' @param font.family the font family for texts
#' @param drop logical to indicate whether all factor levels not used in the data will automatically be dropped. If FALSE (by default), all factor levels will be shown, regardless of whether or not they appear in the data
#' @param details how to label. It can be one of 'name', 'name_FDR' (FDR/OR also appended to the name), and 'name_FDR_members' (FDR/OR plus members appended to the name; in this case, treemap.grow and treemap.reflow is forced to be true)
#' @param caption logical to indicate whether the caption is shown on the bottom-right
#' @param treemap.grow logical to indicate whether text will be grown as well as shrunk to fill the box
#' @param treemap.reflow logical to indicate whether text will be reflowed (wrapped) to better fit the box
#' @param treemap.place where inside the box to place the text. Default is "centre"; other options are "bottom", "top", "topleft", "topright", etc
#' @param treemap.color the color of the text
#' @param treemap.fontface the fontface of the text
#' @param treemap.min.size the minimum font size, in points. If provided, text that would need to be shrunk below this size to fit the box will not be drawn. Defaults to 4 pt
#' @param area which statistics will be used for the area. It can be "adjp" for adjusted p value (FDR) and "or" for odds ratio
#' @param area.fill which statistics will be used for the area fill color. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichTreemap.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' library(treemapify)
#' 
#' # provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # 1) Gene-based enrichment analysis using REACTOME pathways
#' # perform enrichment analysis
#' eTerm <- xEnricherGenes(data, ontology="REACTOME", RData.location=RData.location)
#' ## forest plot of enrichment results
#' gp <- xEnrichTreemap(eTerm, top_num=20, FDR.cutoff=0.05, treemap.reflow=F, treemap.place="topleft")
#' 
#' # 2) Gene-based enrichment analysis using ontologies (REACTOME and GOMF)
#' # perform enrichment analysis
#' ls_eTerm <- xEnricherGenesAdv(data, ontologies=c("REACTOME","GOMF"), RData.location=RData.location)
#' ## forest plot of enrichment results
#' gp <- xEnrichTreemap(ls_eTerm, FDR.cutoff=0.1)
#' }

xEnrichTreemap <- function(eTerm, top_num=10, FDR.cutoff=0.05, CI.one=T, colormap="spectral.top", ncolors=64, zlim=NULL, barwidth=NULL, barheight=0.5, wrap.width=NULL, font.family="sans", drop=F, details=c("name","name_FDR","name_FDR_members"), caption=T, treemap.grow=F, treemap.reflow=F, treemap.place="topleft", treemap.color="black", treemap.fontface="bold.italic", treemap.min.size=4, area=c("adjp","or"), area.fill=c("or","adjp","zscore"))
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    details <- match.arg(details)
    area <- match.arg(area)
    area.fill <- match.arg(area.fill)
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(class(eTerm)=='eTerm'){
		## when 'auto', will keep the significant terms
		df <- xEnrichViewer(eTerm, top_num="all")
		
		############
		if(!CI.one){
			ind <- which(df$CIl>1 | df$CIu<1)
			df <- df[ind,]
		}
		############
		
		if(top_num=='auto'){
			top_num <- sum(df$adjp<FDR.cutoff)
			if(top_num <= 1){
				top_num <- 10
			}
		}
		df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="or", details=T)
		df$group <- 'group'
		df$ontology <- 'ontology'
		
	}else if(class(eTerm)=='ls_eTerm' | class(eTerm)=='data.frame'){
	
		if(class(eTerm)=='ls_eTerm'){
			## when 'auto', will keep the significant terms
			df <- eTerm$df
			
		}else if(class(eTerm)=='data.frame'){
			
			if(all(c('group','ontology','name','zscore','adjp','or','CIl','CIu','nOverlap','members_Overlap') %in% colnames(eTerm))){
				df <- eTerm[,c('group','ontology','name','zscore','adjp','or','CIl','CIu','nOverlap','members_Overlap')]
			
			}else{
				details <- 'name'
				
				if(all(c('group','ontology','name','zscore','adjp','or','CIl','CIu') %in% colnames(eTerm))){
					df <- eTerm[,c('group','ontology','name','zscore','adjp','or','CIl','CIu')]
				
				}else if(all(c('group','name','zscore','adjp','or','CIl','CIu') %in% colnames(eTerm))){
					df <- eTerm[,c('group','name','zscore','adjp','or','CIl','CIu')]
					df$ontology <- 'ontology'
			
				}else if(all(c('ontology','name','zscore','adjp','or','CIl','CIu') %in% colnames(eTerm))){
					df <- eTerm[,c('ontology','name','zscore','adjp','or','CIl','CIu')]
					df$group <- 'group'
			
				}else if(all(c('name','zscore','adjp','or','CIl','CIu') %in% colnames(eTerm))){
					df <- eTerm[,c('name','zscore','adjp','or','CIl','CIu')]
					df$group <- 'group'
					df$ontology <- 'ontology'
			
				}else{
					warnings("The input data.frame does not contain required columns: c('group','ontology','name','zscore','adjp','or','CIl','CIu').\n")
					return(NULL)
				}
			}
			

			
		}
		
		## columns are ordered as indicated by inputs
		df$group <- factor(df$group, levels=unique(df$group))
		
		############
		if(!CI.one){
			ind <- which(df$CIl>1 | df$CIu<1)
			df <- df[ind,]
		}
		############
		
		or <- group <- ontology <- rank <- NULL
		df <- df %>% dplyr::arrange(-or)
		if(top_num=='auto'){
			df <- subset(df, df$adjp<FDR.cutoff)
		}else{
			top_num <- as.integer(top_num)
			df_tmp <- as.data.frame(df %>% dplyr::group_by(group,ontology) %>% dplyr::group_by(rank=rank(-or,decreasing=T),add=TRUE) %>% dplyr::filter(rank<=top_num))
			df <- subset(df, df$name %in% df_tmp$name)
			df <- subset(df, df$adjp<FDR.cutoff)
		}
		
	}

	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}
	
	name <- fdr <- or <- CIl <- CIu <- NULL
	group <- ontology <- NULL
	label <- NULL
	lor <- zscore <- NULL
	
	df$fdr <- -log10(df$adjp)
	df$lor <- log2(df$or)
	
	if(area.fill=='adjp'){
		title <- expression(-log[10]("FDR"))
		if(is.null(zlim)){
			tmp <- df$fdr
			zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
		}
		df$fdr[df$fdr<=zlim[1]] <- zlim[1]
		df$fdr[df$fdr>=zlim[2]] <- zlim[2]
	}else if(area.fill=='or'){
		title <- expression(log[2]("odds ratio"))
		if(is.null(zlim)){
			tmp <- df$lor
			zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
		}
		df$lor[df$lor<=zlim[1]] <- zlim[1]
		df$lor[df$lor>=zlim[2]] <- zlim[2]
	}else if(area.fill=='zscore'){
		title <- "Z-score"
		if(is.null(zlim)){
			tmp <- df$zscore
			zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
		}
		df$zscore[df$zscore<=zlim[1]] <- zlim[1]
		df$zscore[df$zscore>=zlim[2]] <- zlim[2]
	}
	
	## order by 'or', 'adjp'
	df <- df[with(df,order(group, ontology, or, fdr)),]
	df$name <- factor(df$name, levels=unique(df$name))
	
	## label
	if(details=='name'){
		df$label <- df$name
	}else if(details=="name_FDR"){
		df$label <- paste0(df$name, '\n[OR=', df$or, '; FDR=', df$adjp, '; n=', df$nOverlap, ']')
	}else if(details=="name_FDR_members"){
		treemap.grow <- T
		treemap.reflow <- T
		df$label <- paste0(df$name, '\n[OR=', df$or, '; FDR=', df$adjp, '; n=', df$nOverlap, ']\n(', df$members_Overlap,')')
	}
	###########################################
	if(area.fill=='adjp'){
		bp <- ggplot(df, aes(area=lor, fill=fdr, label=label)) 
	}else if(area.fill=='or'){
		bp <- ggplot(df, aes(area=fdr, fill=lor, label=label))
	}else if(area.fill=='zscore'){
		if(area=='adjp'){
			bp <- ggplot(df, aes(area=fdr, fill=zscore, label=label))
		}else if(area=='or'){
			bp <- ggplot(df, aes(area=lor, fill=zscore, label=label))
		}
	}
	bp <- bp + treemapify::geom_treemap() + treemapify::geom_treemap_text(fontface=treemap.fontface, color=treemap.color, place=treemap.place, grow=treemap.grow, reflow=treemap.reflow, min.size=treemap.min.size)
	bp <- bp + theme_bw() + theme(legend.position="bottom")
	bp <- bp + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=title,title.position="left",barwidth=barwidth,barheight=barheight,draw.ulim=FALSE,draw.llim=FALSE))
	
	## caption
    if(caption){
    	if(area=='adjp'){
			bp <- bp + labs(caption="The area is proportional to FDR")
		}else if(area=='or'){
			bp <- bp + labs(caption="The area is proportional to odds ratio")
		}
		bp <- bp + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }

	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	# facet_grid: partitions a plot into a matrix of panels
	## group (columns), ontology (rows)
	ngroup <- length(unique(df$group))
	nonto <- length(unique(df$ontology))
	if(ngroup!=1 | nonto!=1){
		scales <- "free_y"
		space <- "free_y"
		
		if(ngroup==1){
			bp <- bp + facet_grid(ontology~., scales=scales, space=space, drop=drop)
		}else if(nonto==1){
			bp <- bp + facet_grid(.~group, scales=scales, space=space, drop=drop)
		}else{
			bp <- bp + facet_grid(ontology~group, scales=scales, space=space, drop=drop)
		}
		
		## strip
		bp <- bp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=8,face="bold.italic"))
	}
	
	invisible(bp)
}

