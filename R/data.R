#' @name hGOBP.gmt
#' @title Homo sapiens GO biological process (BP) genesets
#' @description BP genesets were downloaded from Gary Bader's
#' lab in July 2022: human_go_bp_no_go_iea_entrezgene.gmt
#' (biological process minus those annoted using evidence code
#' inferred by electronic annotation (IEA))
#' @format gmt file: a list of BP genesets where names = GO term
#' @source \url{http://baderlab.org/GeneSets}
NULL
#' @name yGOBP.gmt
#' @title Saccharomyces GO biological process (BP) genesets
#' @description BP genesets were downloaded in January 2022
#' from the Gene Ontolog Consortium (link below).
#' @format gmt file: a list of BP genesets where names = GO term
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
NULL
#' @name dfGOBP
#' @title Dataframe listing Saccharomyces cerevisiae GOIDs and terms,
#' geneSetSizes, geneSet members
#' @description Dataframe GOID terms correspond to names(yGOBP)
#' annotations were downloaded in January 2022 from Gene Ontolog
#' Consortium (link below)
#' @format Dataframe listing GOIDs and terms
#'      \describe{
#'          \item{GOID}{GO ID}
#'          \item{term}{GO term}
#'          \item{gene}{GO genesets}
#'          \item{geneSetSize}{GO geneset size}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
NULL
#' Chemogenomic Profilimg data from Saccharomyces cerevisiae
#' @name sampleFITdata
#' @title sample fitness data from Saccharomyces cerevisiae
#' @keywords chemogenomics, competitive fitness profiling
#' @description HOP chemogenomic screening data from G. Giaever & C. Nislow labs
#' @format numeric 4956 x 2 matrix, genes in rows,sample in columns
#' @author G.Giaever & C.Nislow
#' @source \url{https://www.ggcnlab.com/}
NULL
