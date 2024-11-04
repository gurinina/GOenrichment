#' Visualize a gene network using visNetwork
#' This function visualizes a gene network from GO enrichment analysis results.
#' It requires pre-processed node and edge data to create an interactive network plot.
#' @param nodes data.frame; contains node information for the network (see `runGORESP` enrichInfo output).
#' @param edges data.frame; contains edge information for the network (see `runGORESP` edgeMat output).
#' @param height Numeric; height of the network plot in pixels. Default is 1300.
#' @param main List; main title settings with `text` and `style`. Default is `list(text = nodes$filename[1], style = "font-size:40px;text-align:center;")`
#' @return No return value, called for side effects. Generates an interactive network visualization in the viewer.
#' @examples
#' # Assuming `nodes` and `edges` are data frames obtained from `runGORESP`
#' runNetwork(nodes, edges, height = 1300, main = list(text = nodes$filename[1], style = "font-size:40px;text-align:center;"))
#' @export
runNetwork <- function(nodes, edges, height = 1300, main = list(text = nodes$filename[1],
    style = "font-size:40px;text-align:center;")) {
    n <- nodes %>% dplyr::arrange(term)
    names <- n$id
    visNetwork(nodes, edges, width = "100%", height = height, main = main) %>%
        visNodes(shadow = list(enabled = TRUE, size = 25)) %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 5, hover = TRUE),
                   nodesIdSelection = list(enabled = TRUE, values = names,
                   style = "width: 400px; height: 31px; font-size: 18px;
                     color: #000066;border: 3px solid #4d88ff;"),
                   selectedBy = list(variable = "FDR", style = "width: 400px;
                     height: 31px; font-size: 18px; color: #000066;
                     border: 3px solid #4d88ff;")) %>%
        visIgraphLayout(type = "full")
}
