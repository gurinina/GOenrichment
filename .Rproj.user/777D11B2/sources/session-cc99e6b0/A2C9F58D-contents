#' Generates a network of GO enrichments.using the output from visSetup
#' @param nodes dataframe returned from visSetup, "nodes" see ?? visSetup
#' @param edges dataframe returned from visSetup, "edges" see ?? visSetup
#' @param height A numeric value
#' @param main character title list, includes font size and color if desired
#' @return network of GO enrichments: The larger a node in the netwark, the more significance the GO enrichment for that term. The thicker the edge connecting twa nodes, the greater the overlap betweem two terms. Bold and enlarged font size indicate the most enriched term for a given cluster.
#' \describe{
#'   \item{querySetFraction}{the fraction of the query set that overlaps with the term set}
#'   \item{geneSetFraction}{the fraction of the term set that overlaps with the query set}
#'   \item{foldEnrichment}{the fold enrichment of the query set with the term genes}
#'   \item{P}{P value estimating the significance with which the query set is enriched with the term genes}
#'   \item{FDR}{FDR value estimating the significance of enrichment }
#'   \item{overlapGenes}{A |-separated list of genes in the overlap of the query set and the term set}
#'   \item{mat}{numeric matrix of fitness data}
#'   \item{coln}{mumeric. the column of the matrix with sample of interest}
#'   \item{curr_exp}{character, name of exp, usually after the sample name so set to colnames(mat)[coln]}
#'   \item{sig}{mumeric, significance threshold of fitness defect score, default is 1}
#'   \item{fdrThresh}{numeric, fdr threshold for significance cutoff of enrichments, default = 0.2}
#'   \item{bp_path}{character, path of gmt file}
#'   \item{bp_input}{gmt file or NULL, **Note: bp_path and bp_input can't both be NULL}
#'   \item{go_path}{character, path of file with GOID and terms path of gmt file alternatively bp_input, gmt file itself,}
#'   \item{go_input}{dataframe with GOID and terms matching bp_input  **Note: go_path and go_input can both be NULL}
#'   \item{minSetSize}{numeric, lower limit on the number of genes in a geneset included in the analysis, default = 5}
#'   \item{maxSetSize}{numeric, upper limit on the number of genes in a geneset included in the analysis,default = 300}
#' @importFrom dplyr arrange desc
#' @examples
#' # Example usage:
#' scoreMat = compSCORE(mat, coln, sig = sig)
#' curr_exp = colnames(mat)[coln]
#' hresp <- GOenrichment::runGORESP(fdrThresh = 0.2, scoreMat = scoreMat,
#'   coln=1, curr_exp = curr_exp, sig = 1,
#'   bp_input = hGOBP.gmt, go_input = NULL, minSetSize = 20,
#'   maxSetSize = 300)
#' @export
\arguments{
\item{nodes}{dataframe returned from visSetup, "nodes" see ?? visSetup}

\item{edges}{dataframe returned from visSetup, "edges" see ?? visSetup}

\item{height}{A numeric value}

\item{main}{character title list, includes font size and color if desired}
}
\value{
network
}
\description{
Uses the output from visSetup to generate a network of gene enrrichments.
}
\details{
The larger a node in the netwark, the more significance the
GO enrichment for that term. The thicker the edge connecting twa nodes,
the greater the overlap betweem twp terms. Bold and enlarged font size
indicate the most enriched term for a given cluster.
}
runNetwork <- function (nodes, edges, height = 1300, main = list(text = nodes$filename[1],
    style = "font-size:40px;text-align:center;"))
{
    n = nodes
    n <- n %>% dplyr::arrange(term)
    names = n$id
    visNetwork(nodes, edges, width = "100%", height = height,
        main = main) %>% visNodes(shadow = list(enabled = T,
        size = 25)) %>% visOptions(highlightNearest = list(enabled = T,
        degree = 5, hover = T), nodesIdSelection = list(enabled = TRUE,
        values = names, style = "width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;"),
        selectedBy = list(variable = "FDR", style = "width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;")) %>%
        visIgraphLayout(type = "full") %>% visEvents(select = "function(nodes) {\nShiny.onInputChange('current_node_id', nodes.nodes);\n;}")
}
