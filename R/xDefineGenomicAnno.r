#' Function to define genomic annotations
#'
#' \code{xDefineGenomicAnno} is supposed to define genomic annotations. It returns an object of class "GenomicRangesList" (GRL).
#'
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in the section 'Note'. Alternatively, the user can also directly provide a customised GR object (or a list of GR objects or a GRL object)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a GRL object
#' @note The genomic annotation data are described below according to the data sources and data types.\cr
#' 1. ENCODE Transcription Factor ChIP-seq data
#' \itemize{
#'  \item{\code{Uniform_TFBS}: a list (690 combinations of cell types and transcription factors) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type per transcription factor.}
#'  \item{\code{ENCODE_TFBS_ClusteredV3}: a list (161 transcription factors) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor, along with a meta-column 'cells' telling cell types associtated with a clustered peak.}
#'  \item{\code{ENCODE_TFBS_ClusteredV3_CellTypes}: a list (91 cell types) of a list (transcription factors) of GenomicRanges objects. Each cell type is a list (transcription factor) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor.}
#' }
#' 2. ENCODE DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{Uniform_DNaseI_HS}: a list (125 cell types) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type.}
#'  \item{\code{ENCODE_DNaseI_ClusteredV3}: an GR object containing clustered peaks, along with a meta-column 'num_cells' telling how many cell types associtated with a clustered peak.}
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
#'  \item{\code{ReMap_Public_TFBS}: a list (1759 combinations of GSE studies and transcription factors and cell types) of GenomicRanges objects; each is an GR object containing identified peaks per GSE study per transcripton factor per cell type.}
#'  \item{\code{ReMap_Encode_TFBS}: a list (1066 combinations of ENCODE transcription factors and cell types) of GenomicRanges objects; each is an GR object containing identified peaks per ENCODE study per transcripton factor per cell type.}
#'  \item{\code{ReMap_PublicAndEncode_TFBS}: a list (2825 combinations of GSE/ENCODE studies and transcription factors and cell types) of GenomicRanges objects; each is an GR object containing identified peaks per GSE/ENCODE study per transcripton factor per cell type.}
#'  \item{\code{ReMap_Public_mergedTFBS}: a list (331 transcription factors under GSE studies) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_Encode_mergedTFBS}: a list (279 transcription factors under ENCODE) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_PublicAndEncode_mergedTFBS}: a list (485 transcription factors under GSE studies and ENCODE) of GenomicRanges objects; each is an GR object containing identified peaks per transcripton factor.}
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
#' 19. GWAS Catalog trait-associated SNPs
#' \itemize{
#'  \item{\code{GWAScatalog_alltraits}: a list (390 traits grouped by EFO) of GenomicRanges objects; each is an GR object containing trait-associated SNPs.}
#'  \item{\code{GWAScatalog_bloodindex}: a list (29 traits grouped by EFO) of GenomicRanges objects; each is an GR object containing trait-associated SNPs.}
#' }
#' @export
#' @seealso \code{\link{xEnrichViewer}}
#' @include xRDataLoader.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' \dontrun{
#' grl <- xDefineGenomicAnno("Uniform_TFBS", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("Uniform_DNaseI_HS", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("FANTOM5_Enhancer_Cell", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("ReMap_Public_TFBS", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("EpigenomeAtlas_15Segments_E029", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("FANTOM5_CAT_Cell", RData.location=RData.location)
#' grl <- xDefineGenomicAnno("GWAScatalog_alltraits", RData.location=RData.location)
#' 
#' # the customised
#' ## a GR object
#' GR.annotation <- grl[[1]]
#' grl_customised <- xDefineGenomicAnno(GR.annotation, RData.location=RData.location)
#' ## a list of GR objects
#' GR.annotation <- lapply(grl[1:2], function(x) x)
#' grl_customised <- xDefineGenomicAnno(GR.annotation, RData.location=RData.location)
#' }

xDefineGenomicAnno <- function(GR.annotation=c(NA,"Uniform_TFBS","ENCODE_TFBS_ClusteredV3","ENCODE_TFBS_ClusteredV3_CellTypes", "Uniform_DNaseI_HS","ENCODE_DNaseI_ClusteredV3","ENCODE_DNaseI_ClusteredV3_CellTypes", "Broad_Histone","SYDH_Histone","UW_Histone","FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_Enhancer_Extensive","FANTOM5_Enhancer","Segment_Combined_Gm12878","Segment_Combined_H1hesc","Segment_Combined_Helas3","Segment_Combined_Hepg2","Segment_Combined_Huvec","Segment_Combined_K562","TFBS_Conserved","TS_miRNA","TCGA", "ReMap_Public_TFBS","ReMap_Encode_TFBS","ReMap_PublicAndEncode_TFBS","ReMap_Public_mergedTFBS","ReMap_Encode_mergedTFBS","ReMap_PublicAndEncode_mergedTFBS","Blueprint_BoneMarrow_Histone","Blueprint_CellLine_Histone","Blueprint_CordBlood_Histone","Blueprint_Thymus_Histone","Blueprint_VenousBlood_Histone","Blueprint_DNaseI","Blueprint_Methylation_hyper","Blueprint_Methylation_hypo","EpigenomeAtlas_15Segments_E029", "EpigenomeAtlas_15Segments_E030", "EpigenomeAtlas_15Segments_E031", "EpigenomeAtlas_15Segments_E032", "EpigenomeAtlas_15Segments_E033", "EpigenomeAtlas_15Segments_E034", "EpigenomeAtlas_15Segments_E035", "EpigenomeAtlas_15Segments_E036", "EpigenomeAtlas_15Segments_E037", "EpigenomeAtlas_15Segments_E038", "EpigenomeAtlas_15Segments_E039", "EpigenomeAtlas_15Segments_E040", "EpigenomeAtlas_15Segments_E041", "EpigenomeAtlas_15Segments_E042", "EpigenomeAtlas_15Segments_E043", "EpigenomeAtlas_15Segments_E044", "EpigenomeAtlas_15Segments_E045", "EpigenomeAtlas_15Segments_E046", "EpigenomeAtlas_15Segments_E047", "EpigenomeAtlas_15Segments_E048", "EpigenomeAtlas_15Segments_E050", "EpigenomeAtlas_15Segments_E051", "EpigenomeAtlas_15Segments_E062", "CpG_anno","Genic_anno", "FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","FANTOM5_CAT_DO","FANTOM5_CAT_EFO","FANTOM5_CAT_HPO","FANTOM5_CAT_MESH","FANTOM5_CAT_PICS"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    
    grl <- NULL
    
    if(class(GR.annotation) == "GenomicRangesList"){
    	grl <- GR.annotation
    	
    }else{
    
		lgr <- NULL
		if(class(GR.annotation) == "list"){
			lgr <- GR.annotation
		
			if(!all(sapply(lgr, function(x) class(x) == "GRanges"))){
				lgr <- NULL
			}
		
		}else if(class(GR.annotation) == "GRanges"){
			lgr <- list(Customised=GR.annotation)
	
		}else{
	
			if(length(GR.annotation)>1){
				GR.annotation <- GR.annotation[1]
			}
	
			if(!is.na(GR.annotation)){
				lgr <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
			}
		}
	
		if(class(lgr) == "list"){
			## Remove null elements in a list
			lgr <- base::Filter(base::Negate(is.null), lgr)
			grl <- GenomicRanges::GRangesList(lgr)	
		}
	
	}
	
	invisible(grl)
}
