#' BP genesets from Homo sapiens
#'
#' GO biological process (BP) genesets from
#' Downloaded from Gary Bader's lab in July 2022
#'
#' @format .gmt file: a list of BP genesets where names = GO term
#'   \describe{
#'      \length{16746}
#' }
#' @source \url{http://baderlab.org/GeneSets}
"hGOBP.gmt"

#' BP genesets from Saccharomyces cerevisiae
#'
#' GO biological process (BP) genesets from
#' Saccharomyces cerevisiae downloaded in January 2022
#' @format .gmt file: a list of BP genesets where names = GO term
#'   \describe{
#'      \length{3227}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"yGOBP.gmt"

#' BP genesets from Saccharomyces cerevisiae
#'
#' GO biological process (BP) genesets from
#' Saccharomyces cerevisiae downloaded in January 2022
#'
#' @format Dataframe listing GOIDs and terms
#'   \describe{
#'      \item{GOID}{GO ID}
#'      \item{term}{GO term}
#'      \item{gene}{GO genesets}
#'      \item{geneSetSize}{GO geneset size}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"dfGOBP"

#' Fitness data from Saccharomyces cerevisiae
#'
#' Experimental data from G. Giaever & C. Nislow labs
#'
#'
#' @format numeric matrix with two columns; 4956 gene names in rows
#'   \describe{
#'      \dim{4956 2}
#' }
#'@source \url{https://www.ggcnlab.com/}
"sampleFITdata"

