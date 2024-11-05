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
visSetup = function(enrichInfo, edgeMat, fontsize = 22, fontface = "Arial") {
  library(dplyr)
  library(igraph)
  library(visNetwork)

  message("dplyr loaded: ", "dplyr" %in% loadedNamespaces())

  n = enrichInfo
  e = edgeMat

  w = which(names(n) == "id")
  coln = (1:ncol(n))[-w]
  n = n[, c(w, coln)]

  if (is.null(e) & !is.null(n)) {
    gr = make_empty_graph(nrow(enrichInfo))
    v = gr
    V(v)$color.background = n$cluster
    v = set_vertex_attr(v, "label", value = n$formattedLabel)

  }


  if (!is.null(e) & !is.null(n)) {
    w = which(names(e) == "label")
    let = graph_from_data_frame(e[, -w], vertices = n, directed = F)
    v = set_vertex_attr(let, "label", value = n$formattedLabel)
    V(v)$color.background = n$cluster
  }

  vis = toVisNetworkData(v)
  vis$nodes = data.frame(vis$nodes, stringsAsFactors = F)
  if (!is.null(vis$edges))
    vis$edges = data.frame(vis$edges, stringsAsFactors = F)
  m = match(vis$nodes$id, n$id)
  vis$nodes$label = n$formattedLabel[m]
  vis$nodes$FDR = n$FDR[m]
  vis$nodes$FDR = signif(vis$nodes$FDR, 5)


  w = which(duplicated(vis$nodes$FDR))
  if (length(w) > 0) {

    vis$nodes <- dplyr::group_by(vis$nodes, FDR)
    vis$nodes <- dplyr::mutate(vis$nodes, jitter = if (dplyr::n() > 1) abs(jitter(FDR)) else FDR)
    vis$nodes <- dplyr::ungroup(vis$nodes)


    w = which(names(vis$nodes) == "FDR")
    vis$nodes = vis$nodes[, -w]
    w = which(names(vis$nodes) == "jitter")
    names(vis$nodes)[w] = "FDR"
    vis$nodes$FDR = signif(vis$nodes$FDR, 6)
  }

  w = which(duplicated(vis$nodes$FDR))



  vis$nodes$term = n$term[m]
  vis$nodes$interSect = n$interSect[m]
  vis$nodes$nQuery= n$nQuery[m]
  vis$nodes$nGenes = n$nGenes[m]
  vis$nodes$geneSetFraction = n$geneSetFraction[m]
  vis$nodes$querySetFraction = n$querySetFraction[m]
  vis$nodes$filename =  n$filename[m]

  vis$nodes$formattedLabel = n$formattedLabel[m]
  vis$nodes$overlapGenes = n$overlapGenes[m]
  vis$nodes$label = vis$nodes$formattedLabel
  vis$nodes$color.border = "black"
  vis$nodes = vis$nodes %>% dplyr::arrange(label)
  if (nrow(vis$edges) > 0)
    vis$edges$color = "black"

  vis$nodes$borderWidthSelected	= 4
  w = which(names(vis$nodes) %in% c("fomattedLabel", "color", "cluster"))
  if (length(w) > 0)  vis$nodes = vis$nodes[, -w]
  vis$nodes$color.highlight.border = "#000066"
  vis$nodes$color.highlight.background = "#c0b3ff"
  vis$nodes$color.hover.background = "#000066"
  vis$nodes$color.hover.border = "#c0b3ff"
  vis$nodes$font.face = fontface
  vis$nodes$shape = "dot"
  vis$nodes$font.size = fontsize
  vis$nodes$font.bold = F
  vis$nodes$borderWidth = 2
  #vis$nodes$vadjust = "mono"
  vis$nodes$borderWidthSelected = 4
  vis$nodes$labelHighlightBold = T
  w = which.min(vis$nodes$FDR)
  if (length(w) > 0) {
    vis$nodes$font.size[w] = fontsize + 4
    vis$nodes$borderWidth[w] = 4
    vis$nodes$font.bold[w] = T
  }
  vis$nodes = vis$nodes %>% dplyr::arrange(FDR)
  vis

}
