#' Function to calculate per base scores given a list of genomic regions in terms of overlaps with genomic annotations
#'
#' \code{xGScoreAdv} is supposed to calculate per base scores for an input list of genomic regions (genome build 19), using genomic annotations (eg genomic segments, active chromatin, transcription factor binding sites/motifs, conserved sites). The per base scores are calculated for overlaps with each genomic annotation. Scores for genomic regions/variants can be constraint/conservation or impact/pathogenicity.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GS.annotation which genomic scores (GS) annotaions used. It can be 'fitCons' (the probability of fitness consequences for point mutations; \url{http://www.ncbi.nlm.nih.gov/pubmed/25599402}), 'phastCons' (the probability that each nucleotide belongs to a conserved element/negative selection [0,1]), 'phyloP' (conservation at individual sites representing -log p-values under a null hypothesis of neutral evolution, positive scores for conservation and negative scores for acceleration), 'mcap' (eliminating a majority of variants with uncertain significance in clinical exomes at high sensitivity: \url{http://www.ncbi.nlm.nih.gov/pubmed/27776117}), and 'cadd' (combined annotation dependent depletion for estimating relative levels of pathogenicity of potential human variants: \url{http://www.ncbi.nlm.nih.gov/pubmed/24487276})
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in the section 'Note'. Beyond pre-built annotation data, the user can specify the customised input. To do so, first save your RData file (a list of GR objects, each is an GR object correponding to an annotation) into your local computer. Then, tell "GR.annotation" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param details logical to indicate whether the detailed information (ie ratio) is returned. By default, it sets to false for no inclusion
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 6 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{o_nBase}: the number of bases overlapped between input regions and annotation regions}
#'  \item{\code{o_GS}: the per base genomic scores for overlaps between input regions and annotation regions}
#'  \item{\code{a_nBase}: the number of bases covered by that annotation; optional, it is only appended when "details" is true}
#'  \item{\code{a_GS}: the per base genomic scores for that annotation; optional, it is only appended when "details" is true}
#'  \item{\code{ratio}: ratio of o_GS divided by a_GS; optional, it is only appended when "details" is true}
#' }
#' @note The genomic annotation data are described below according to the data sources and data types.\cr
#' 1. ENCODE Transcription Factor ChIP-seq data
#' \itemize{
#'  \item{\code{Uniform_TFBS}: a list (690 combinations of cell types and transcription factors) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type per transcription factor.}
#'  \item{\code{ENCODE_TFBS_ClusteredV3_CellTypes}: a list (91 cell types) of a list (transcription factors) of GenomicRanges objects. Each cell type is a list (transcription factor) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor.}
#' }
#' 2. ENCODE DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{Uniform_DNaseI_HS}: a list (125 cell types) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type.}
#'  \item{\code{ENCODE_DNaseI_ClusteredV3_CellTypes}: a list (125 cell types) of GenomicRanges objects; each is an GR object containing clustered peaks per cell type.}
#' }
#' 3. ENCODE Histone Modification ChIP-seq data from different sources
#' \itemize{
#'  \item{\code{Broad_Histone}: a list (156 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/Broad Institute.}
#'  \item{\code{SYDH_Histone}: a list (29 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/Stanford/Yale/Davis/Harvard.}
#'  \item{\code{UW_Histone}: a list (172 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/University of Washington.}
#' }
#' 4. FANTOM5 expressed enhancer atlas
#' \itemize{
#'  \item{\code{FANTOM5_Enhancer_Cell}: a list (71 cell types) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_Enhancer_Tissue}: a list (41 tissues) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a tissue.}
#'  \item{\code{FANTOM5_Enhancer_Extensive}: a list (5 categories of extensitive enhancers) of GenomicRanges objects; each is an GR object containing extensitive enhancers. They are: "Extensive_ubiquitous_enhancers_cells" for ubiquitous enhancers expressed over the entire set of cell types; "Extensive_ubiquitous_enhancers_organs" for ubiquitous enhancers expressed over the entire set of tissues; "Extensive_enhancers_tss_associations" for TSS-enhancer associations(RefSeq promoters only); "Extensive_permissive_enhancers" and "Extensive_robust_enhancers" for permissive and robust enhancer sets.}
#'  \item{\code{FANTOM5_Enhancer}: a list (117 cell types/tissues/categories) of GenomicRanges objects; each is an GR object.}
#' }
#' 5. ENCODE combined (ChromHMM and Segway) Genome Segmentation data
#' \itemize{
#'  \item{\code{Segment_Combined_Gm12878}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line GM12878 (a lymphoblastoid cell line).}
#'  \item{\code{Segment_Combined_H1hesc}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line H1-hESC (H1 human embryonic stem cells).}
#'  \item{\code{Segment_Combined_Helas3}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HeLa S3.}
#'  \item{\code{Segment_Combined_Hepg2}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HepG2 (liver hepatocellular carcinoma).}
#'  \item{\code{Segment_Combined_Huvec}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HUVEC (Human Umbilical Vein Endothelial Cells).}
#'  \item{\code{Segment_Combined_K562}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line K562 (human erythromyeloblastoid leukemia cell line).}
#' }
#' 6. Conserved TFBS
#' \itemize{
#'  \item{\code{TFBS_Conserved}: a list (245 PWM) of GenomicRanges objects; each is an GR object containing human/mouse/rat conserved TFBS for each PWM.}
#' }
#' 7. TargetScan miRNA regulatory sites
#' \itemize{
#'  \item{\code{TS_miRNA}: a list (153 miRNA) of GenomicRanges objects; each is an GR object containing miRNA regulatory sites for each miRNA.}
#' }
#' 8. TCGA exome mutation data
#' \itemize{
#'  \item{\code{TCGA}: a list (11 tumor types) of GenomicRanges objects; each is an GR object containing exome mutation across tumor patients of the same tumor type.}
#' }
#' 9. ReMap integration of transcription factor ChIP-seq data (publicly available and ENCODE)
#' \itemize{
#'  \item{\code{ReMap_Public_TFBS}: a list (395 combinations of GSE studies and transcription factors and cell types) of GenomicRanges objects; each is an GR object containing identified peaks per GSE study per transcripton factor per cell type.}
#'  \item{\code{ReMap_Public_mergedTFBS}: a list (131 transcription factors under GSE studies) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_PublicAndEncode_mergedTFBS}: a list (237 transcription factors under GSE studies and ENCODE) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_Encode_TFBS}: a list (155 transcription factors under ENCODE) of GenomicRanges objects; each is an GR object containing identified peaks per transcripton factor.}
#' }
#' 10. Blueprint Histone Modification ChIP-seq data
#' \itemize{
#'  \item{\code{Blueprint_BoneMarrow_Histone}: a list (132 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from bone marrow).}
#'  \item{\code{Blueprint_CellLine_Histone}: a list (38 combinations of histone modifications and cell lines) of GenomicRanges objects; each is an GR object containing identified peaks per histone per cell line.}
#'  \item{\code{Blueprint_CordBlood_Histone}: a list (126 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from cord blood).}
#'  \item{\code{Blueprint_Thymus_Histone}: a list (5 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from thymus).}
#'  \item{\code{Blueprint_VenousBlood_Histone}: a list (296 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from venous blood).}
#' }
#' 11. BLUEPRINT DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{Blueprint_DNaseI}: a list (36 samples) of GenomicRanges objects; each is an GR object containing identified peaks per sample.}
#' }
#' 12. BLUEPRINT DNA Methylation data
#' \itemize{
#'  \item{\code{Blueprint_Methylation_hyper}: a list (206 samples) of GenomicRanges objects; each is an GR object containing hyper-methylated CpG regions per sample.}
#'  \item{\code{Blueprint_Methylation_hypo}: a list (206 samples) of GenomicRanges objects; each is an GR object containing hypo-methylated CpG regions per sample.}
#' }
#' 13. Roadmap Epigenomics Core 15-state Genome Segmentation data for primary cells (blood and T cells)
#' \itemize{
#' \item{\code{EpigenomeAtlas_15Segments_E033}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E033 (Primary T cells from cord blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E034}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E034 (Primary T cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E037}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E037 (Primary T helper memory cells from peripheral blood 2).}
#' \item{\code{EpigenomeAtlas_15Segments_E038}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E038 (Primary T helper naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E039}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E039 (Primary T helper naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E040}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E040 (Primary T helper memory cells from peripheral blood 1).}
#' \item{\code{EpigenomeAtlas_15Segments_E041}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E041 (Primary T helper cells PMA-I stimulated).}
#' \item{\code{EpigenomeAtlas_15Segments_E042}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E042 (Primary T helper 17 cells PMA-I stimulated).}
#' \item{\code{EpigenomeAtlas_15Segments_E043}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E043 (Primary T helper cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E044}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E044 (Primary T regulatory cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E045}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E045 (Primary T cells effector/memory enriched from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E047}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E047 (Primary T killer naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E048}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E048 (Primary T killer memory cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E062}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E062 (Primary mononuclear cells from peripheral blood).}
#' }
#' 14. Roadmap Epigenomics Core 15-state Genome Segmentation data for primary cells (HSC and B cells)
#' \itemize{
#' \item{\code{EpigenomeAtlas_15Segments_E029}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E029 (Primary monocytes from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E030}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E030 (Primary neutrophils from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E031}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E031 (Primary B cells from cord blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E032}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E032 (Primary B cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E035}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E035 (Primary hematopoietic stem cells).}
#' \item{\code{EpigenomeAtlas_15Segments_E036}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E036 (Primary hematopoietic stem cells short term culture).}
#' \item{\code{EpigenomeAtlas_15Segments_E046}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E046 (Primary Natural Killer cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E050}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E050 (Primary hematopoietic stem cells G-CSF-mobilized Female).}
#' \item{\code{EpigenomeAtlas_15Segments_E051}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E051 (Primary hematopoietic stem cells G-CSF-mobilized Male).}
#' }
#' 15. CpG annotation
#' \itemize{
#'  \item{\code{CpG_anno}: a list (4 categories) of GenomicRanges objects; each is an GR object. They are exclusive, including (in order) "CpG_islands", "CpG_shores" (2Kb upstream/downstream from the ends of the CpG islands), "CpG_shelves" (2Kb upstream/downstream of the farthest upstream/downstream limits of the CpG shores), and "CpG_inter" (the remaining inter-CGI genomic regions 'open sea'). }
#' }
#' 16. Genic annotation
#' \itemize{
#'  \item{\code{Genic_anno}: a list (12 categories) of GenomicRanges objects; each is an GR object. They are not exclusively, including "Genic_1to5kb" (1-5Kb upstream of TSS), "Genic_promoters" (1Kb upstream of TSS), "Genic_5UTRs", "Genic_firstexons" (first exons), "Genic_exons", "Genic_exonintronboundaries", "Genic_introns", "Genic_intronexonboundaries", "Genic_cds", "Genic_3UTRs", "Genic_intergenic" (the intergenic regions exclude the previous list of annotations), and "Genic_lncRNA" (GENCODE long non-coding RNA (lncRNA) transcripts). }
#' }
#' 17. FANTOM5 sample-ontology-enriched CAT genes
#' \itemize{
#'  \item{\code{FANTOM5_CAT_Cell}: a list (173 cell types) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_CAT_Tissue}: a list (174 tissues) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a tissue.}
#' }
#' 18. FANTOM5 trait-associated CAT genes
#' \itemize{
#'  \item{\code{FANTOM5_CAT_DO}: a list (299 traits grouped by disease ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_EFO}: a list (93 traits grouped by experiment factor ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_HPO}: a list (176 traits grouped by human phenotype ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_MESH}: a list (210 traits grouped by Medical Subject Headings) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_PICS}: a list (39 traits grouped by PICS dieases) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#' }
#' @export
#' @seealso \code{\link{xGScore}}
#' @include xGScoreAdv.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS
#' data <- ImmunoBase$AS$variant
#'
#' # b) in terms of overlaps with genomic segments (Primary monocytes from peripheral blood)
#' ## fitness consequence score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="fitCons", GR.annotation="EpigenomeAtlas_15Segments_E029", RData.location=RData.location)
#' ## phastCons conservation score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="phastCons", GR.annotation="EpigenomeAtlas_15Segments_E029", RData.location=RData.location)
#' 
#' # c) in terms of overlaps with genic annotations
#' ## phyloP conservation score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="phyloP", GR.annotation="Genic_anno", RData.location=RData.location)
#' }

xGScoreAdv <- function(data, format=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GS.annotation=c("fitCons","phastCons","phyloP","mcap","cadd"), GR.annotation=c(NA,"Uniform_TFBS","ENCODE_TFBS_ClusteredV3","ENCODE_TFBS_ClusteredV3_CellTypes", "Uniform_DNaseI_HS","ENCODE_DNaseI_ClusteredV3","ENCODE_DNaseI_ClusteredV3_CellTypes", "Broad_Histone","SYDH_Histone","UW_Histone","FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_Enhancer_Extensive","FANTOM5_Enhancer","Segment_Combined_Gm12878","Segment_Combined_H1hesc","Segment_Combined_Helas3","Segment_Combined_Hepg2","Segment_Combined_Huvec","Segment_Combined_K562","TFBS_Conserved","TS_miRNA","TCGA", "ReMap_Public_TFBS","ReMap_Public_mergedTFBS","ReMap_PublicAndEncode_mergedTFBS","ReMap_Encode_TFBS", "Blueprint_BoneMarrow_Histone","Blueprint_CellLine_Histone","Blueprint_CordBlood_Histone","Blueprint_Thymus_Histone","Blueprint_VenousBlood_Histone","Blueprint_DNaseI","Blueprint_Methylation_hyper","Blueprint_Methylation_hypo","EpigenomeAtlas_15Segments_E029", "EpigenomeAtlas_15Segments_E030", "EpigenomeAtlas_15Segments_E031", "EpigenomeAtlas_15Segments_E032", "EpigenomeAtlas_15Segments_E033", "EpigenomeAtlas_15Segments_E034", "EpigenomeAtlas_15Segments_E035", "EpigenomeAtlas_15Segments_E036", "EpigenomeAtlas_15Segments_E037", "EpigenomeAtlas_15Segments_E038", "EpigenomeAtlas_15Segments_E039", "EpigenomeAtlas_15Segments_E040", "EpigenomeAtlas_15Segments_E041", "EpigenomeAtlas_15Segments_E042", "EpigenomeAtlas_15Segments_E043", "EpigenomeAtlas_15Segments_E044", "EpigenomeAtlas_15Segments_E045", "EpigenomeAtlas_15Segments_E046", "EpigenomeAtlas_15Segments_E047", "EpigenomeAtlas_15Segments_E048", "EpigenomeAtlas_15Segments_E050", "EpigenomeAtlas_15Segments_E051", "EpigenomeAtlas_15Segments_E062", "CpG_anno","Genic_anno", "FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","FANTOM5_CAT_DO","FANTOM5_CAT_EFO","FANTOM5_CAT_HPO","FANTOM5_CAT_MESH","FANTOM5_CAT_PICS"), details=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata_dev")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    GS.annotation <- match.arg(GS.annotation)
    
	if(class(GR.annotation) == "GRanges"){
		###################################
		## now GR.annotation can be directly provided as a GR object
		###################################
		aGR <- GR.annotation
	}else{
		if(is.na(GR.annotation)){
			return(NULL)
		}else{
			if(length(GR.annotation)>1){
				message("\tONLY the first specified annotation RData will be used!\n")
				GR.annotation <- GR.annotation[1]
			}
			aGR <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
			if(is.null(aGR)){
				stop("Your specified annotation RData does not exist!\n")
			}
		}
	}
    
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, maxgap=-1L, minoverlap=0L){
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits <- qGR[hits[,1]]
		shits <- sGR[hits[,2]]

		oGR <- IRanges::pintersect(qhits, shits, ignore.strand=T)
		#IRanges::reduce(oGR)
	}
	#####################################
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	ls_df <- lapply(1:length(aGR), function(i){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Calculating '%s' for the annotation '%s' (%s) ...", GS.annotation, names(aGR)[i], as.character(now)), appendLF=TRUE)
		}
		gr <- aGR[[i]]
		
		## overlaps
		gr_overlap <- mergeOverlaps(dGR, gr)
		gr_score_overlap <- xGScore(data=gr_overlap, format="GRanges", GS.annotation=GS.annotation, scoring.scheme="sum", RData.location=RData.location, verbose=FALSE)
		ind <- which(!is.na(gr_score_overlap$GScore))
		o_nBase <- sum(as.numeric(IRanges::width(gr_score_overlap[ind])))
		o_GScore <- sum(gr_score_overlap$GScore[ind])
		o_GScore_per_base <- o_GScore / o_nBase
		#o_GScore_per_base[is.na(o_GScore_per_base)] <- 0
		
		if(details){
			## annotation
			gr_score_bg <- xGScore(data=gr, format="GRanges", GS.annotation=GS.annotation, scoring.scheme="sum", RData.location=RData.location, verbose=FALSE)
			ind <- which(!is.na(gr_score_bg$GScore))
			a_nBase <- sum(as.numeric(IRanges::width(gr_score_bg[ind])))
			a_GScore <- sum(gr_score_bg$GScore[ind])
			a_GScore_per_base <- a_GScore / a_nBase
			
			## output
			res <- data.frame(name=names(aGR)[i], o_nBase, o_GS=o_GScore_per_base, a_nBase, a_GS=a_GScore_per_base, ratio=o_GScore_per_base/a_GScore_per_base, stringsAsFactors=F)
			
		}else{
			## output
			res <- data.frame(name=names(aGR)[i], o_nBase, o_GS=o_GScore_per_base, stringsAsFactors=F)
		}
		
		return(res)
	})
	res_df <- do.call(rbind, ls_df)
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
	
	invisible(res_df)
}
