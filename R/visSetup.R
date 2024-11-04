#' Prepare network data for gene visualization using visNetwork
#' This function processes the output of `runGORESP` (`enrichInfo` and `edgeMat`) to create formatted nodes and edges for use in visNetwork.
#' @param enrichInfo Data frame of enrichment information, as produced by `runGORESP`.
#' @param edgeMat Data frame of edges, showing relationships between nodes, as produced by `runGORESP`.
#' @param fontsize Numeric; font size for node labels. Default is 22.
#' @param fontface Character; font face for node labels. Default is "Arial".
#' @return A list with two formatted data frames: `nodes` (containing attributes for each node) and `edges` (representing relationships between nodes), both suitable for `visNetwork` plotting.
#' @examples
#' visSetup(enrichInfo, edgeMat)
#' @export
visSetup <- function(enrichInfo, edgeMat, fontsize = 22, fontface = "Arial") {
    n <- enrichInfo
    e <- edgeMat
    w <- which(names(n) == "id")
    coln <- (1:ncol(n))[-w]
    n <- n[, c(w, coln)]
    if (is.null(e) & !is.null(n)) {
        gr <- igraph::make_empty_graph(nrow(enrichInfo))
        v <- gr
        igraph::V(v)$color.background <- n$cluster
        v <- igraph::set_vertex_attr(v, "label", value = n$formattedLabel)
    }
    if (!is.null(e) & !is.null(n)) {
        w <- which(names(e) == "label")
        let <- igraph::graph_from_data_frame(e[, -w], vertices = n, directed = FALSE)
        v <- igraph::set_vertex_attr(let, "label", value = n$formattedLabel)
        igraph::V(v)$color.background <- n$cluster
    }
    vis <- visNetwork::toVisNetworkData(v)
    vis$nodes <- data.frame(vis$nodes, stringsAsFactors = FALSE)
    if (!is.null(vis$edges)) {
        vis$edges <- data.frame(vis$edges, stringsAsFactors = FALSE)
    }

    # Map enrichment information to visNetwork nodes
    m <- match(vis$nodes$id, n$id)
    vis$nodes$label <- n$formattedLabel[m]
    vis$nodes$FDR <- signif(n$FDR[m], 5)
    vis$nodes <- dplyr::group_by(vis$nodes, FDR)
    vis$nodes <- dplyr::mutate(vis$nodes, jitter = if_else(n() > 1, abs(jitter(FDR)), FDR))
    vis$nodes <- dplyr::ungroup(vis$nodes)
    vis$nodes$label <- vis$nodes$formattedLabel
    vis$nodes$color.border <- "black"
    vis$nodes <- dplyr::arrange(vis$nodes, label)

    # Edge formatting
    if (nrow(vis$edges) > 0) vis$edges$color <- "black"
    vis$nodes$borderWidthSelected <- 4
    vis$nodes$color.highlight.border <- "#000066"
    vis$nodes$color.highlight.background <- "#c0b3ff"
    vis$nodes$color.hover.background <- "#000066"
    vis$nodes$color.hover.border <- "#c0b3ff"
    vis$nodes$font.face <- fontface
    vis$nodes$shape <- "dot"
    vis$nodes$font.size <- fontsize
    vis$nodes$font.bold <- FALSE
    vis$nodes$borderWidth <- 2
    vis$nodes$borderWidthSelected <- 4
    vis$nodes$labelHighlightBold <- TRUE

    # Emphasize the node with the lowest FDR
    w <- which.min(vis$nodes$FDR)
    if (length(w) > 0) {
        vis$nodes$font.size[w] <- fontsize + 4
        vis$nodes$borderWidth[w] <- 4
        vis$nodes$font.bold[w] <- TRUE
    }

    vis$nodes <- dplyr::arrange(vis$nodes, FDR)
    vis
}
