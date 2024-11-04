
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
