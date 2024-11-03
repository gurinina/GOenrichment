#' GO enrichment analysis using the hypergeometric test
#' This function performs a GO enrichment analysis using the hypergeometric test for the set of query genes passing the user-defined fitness score threshold compared to the genes in the background set (universe).
#' The hypergeometric test requires both a list of selected genes (e.g. your DGE genes) and a “universe” list (e.g. all genes represented that were tested for differential expression), all represented by their “SYMBOL” gene ID.
#' @param scoreMat A data frame, see ?runGORESP enrichInfo output for explanation of these columns
#' @param coln A data frame, see ?runGORESP edgeMat output for explanation of these columns
#' @param curr_exp fontsize for node labels; default = 22
#' @param sig numeric, significance threshold of fitness defect score, default is 1
#' @param fdrThresh numeric, fdr threshold for significance cutoff of enrichments, default = 0.2
#' @param bp_path character, path of gmt file
#' @param bp_input gmt file or NULL, **Note: bp_path and bp_input can't both be NULL
#' @param go_path character, path of file with GOID and terms path of gmt file;
#' alternatively bp_input, gmt file itself,
#' @param go_input dataframe with GOID and terms matching bp_input **Note: go_path and go_input can both be NULL
#' @param minSetSize numeric, lower limit on the number of genes in a geneset included in the analysis,
#' default = 5
#' @param maxSetSize numeric, upper limit on the number of genes in a geneset included in the analysis,
#' default = 300
#' @param uniSize The number of genes in the universe or length(uni)
#' @return RETURNS two dataframes: enrichInfo with enrichment results and edgeMat with information about
#' geneset term overlap and clusters.
#' dataframe includes the columns: filename,term, nGenes, nQuery, nOverlap,
#' columns enrichInfo:
#' filename
#' GOID: GO ID
#' term: GO term
#' nGenes: # genes in the geneset nQuery: # genes in the query
#' nOverlap: # genes overlapping the query and the geneset
#' querySetFraction: the fraction of the query set (the set that passes the significance threshold) that overlaps with the term set
#' geneSetFraction: the fraction of the term set that overlaps with the query set
#' foldEnrichment = the fold enrichment of the query set over the
#' bgRate, where bgRate is nGenes in geneSet/genes in the universe
#' P = P value estimating the significance with which the query set is enriched with the term genes
#' FDR = FDR value estimating the significance of enrichment
#' overlapGenes = a |-separated list of genes that overlap the query set and the geneSet
#' if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#' maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
#' cluster: cluster color based on overlap of term ID, set as 50
#' id: cluster ID size: node size for downstream plotting based on FDR score formatted
#' label: for downstream plotting
#' columns edgeMat:
#' source id from enrichInfo
#' target id from enrichInfo
#' overlapCoeff is defined as a function of all possible pairs of geneSets
#' overlapCoeff = function(gsPairList) {
#' length(intersect(gsPairList[[1]], gsPairList[[2]]))/min(length(gsPairList[[1]]),
#' length(gsPairList[[2]]))}
#' width
#' label These map GO terms from the enrichInfo data.frame serving to show the overlap between terms for downstream visualization purposes and for understanding the relationship between terms
#' @return A data frame with three columns:
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
runGORESP <- function (scoreMat, coln, curr_exp = colnames(mat)[coln], sig = 1,
    fdrThresh = 0.2, bp_path = NULL, bp_input = NULL, go_path = NULL,
    go_input = NULL, minSetSize = 5, maxSetSize = 300){

    CLUST.COL <- c("#FF00CC", "#33CCFF", "#33CC00", "#9900FF",
        "#FF9900", "#FFFF00", "#FFCCFF", "#FF0000", "#006600",
        "#009999", "#CCCC00", "#993300", "#CC99CC", "#6699CC",
        "#CCCCFF", "#FFCC99", "#9966FF", "#CC6600", "#CCFFFF",
        "#99CC00", "#FF99FF", "#0066FF", "#66FFCC", "#99CCFF",
        "#9999CC", "#CC9900", "#CC33FF", "#006699", "#F5DF16",
        "#B5185E", "#99FF00", "#00FFFF", "#990000", "#CC0000",
        "#33CCCC", "#CC6666", "#996600", "#9999FF", "#3366FF")
    prunedCol <- "#BEBEBE"
    getUniquePairs = function(maxVal) {
        firstI <- rep(1:(maxVal - 1), (maxVal - 1):1)
        secondI <- sapply(2:maxVal, function(x) {
            x:maxVal
        })
        cbind(firstI, unlist(secondI))
    }
    hyperG = function(querySet, geneSets, uni, scoreMat, minSetSize = 5,
        maxSetSize = 300, uniSize = NA) {
        if (!is.null(uni)) {
            geneSets <- lapply(geneSets, intersect, uni)
            lens <- sapply(geneSets, length)
            geneSets <- geneSets[lens >= minSetSize & lens <=
                maxSetSize]
            uniSize <- length(uni)
        }
        if (!is.null(scoreMat)) {
            scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),
                ]
            if (!is.null(uni)) {
                i <- match(uni, scoreMat$gene)
                scoreMat <- scoreMat[sort(i[!is.na(i)]), ]
            }
            scoreMat$score <- round(scoreMat$score, 2)
        }
        enrichInfo <- sapply(geneSets, function(geneSet) {
            overlapSet <- intersect(querySet, geneSet)
            pVal <- stats::phyper(length(overlapSet) - 1, length(geneSet),
                uniSize - length(geneSet), length(querySet),
                lower.tail = F)
            if (length(overlapSet) > 0) {
                overlapSet <- sort(overlapSet)
            }
            overlapSize <- length(overlapSet)
            if (is.null(scoreMat)) {
                maxScore <- NA
            }
            else {
                i <- sort(match(overlapSet, scoreMat$gene))
                maxScore <- scoreMat$score[i[1]]
                overlapSet <- paste(scoreMat$gene[i], "(", scoreMat$score[i],
                  ")", sep = "")
            }
            overlapSet <- paste(overlapSet, collapse = "|")
            bgRate <- length(geneSet)/uniSize
            foldEnrich <- overlapSize/length(querySet)/bgRate
            c(overlapSet, overlapSize/length(geneSet), foldEnrich,
                pVal, maxScore, overlapSize/length(querySet))
        })
        enrichInfo <- t(enrichInfo)
        enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[,
            6]), geneSetFraction = as.numeric(enrichInfo[, 2]),
            foldEnrichment = as.numeric(enrichInfo[, 3]), P = as.numeric(enrichInfo[,
                4]), FDR = stats::p.adjust(as.numeric(enrichInfo[,
                4]), method = "BH"), overlapGenes = enrichInfo[,
                1], maxOverlapGeneScore = as.numeric(enrichInfo[,
                5]), stringsAsFactors = F)
        rownames(enrichCol) <- NULL
        enrichCol = enrichCol[order(enrichCol$FDR), ]
    }
    overlapCoeff = function(gsPairList) {
        length(intersect(gsPairList[[1]], gsPairList[[2]]))/min(length(gsPairList[[1]]),
            length(gsPairList[[2]]))
    }
    clusterEnrich = function(enrichInfo, geneSets, fdrThresh = 0.1,
        overlapThresh = 0.5, goTable = NULL) {
        nodeSizeRange <- c(10, 40)
        prunedCol <- "#BEBEBE"
        labelWidth <- 20
        edgeWidthRange <- c(1, 5)
        overlapCoeffRange <- c(overlapThresh, 1)
        enrichInfo <- enrichInfo[enrichInfo$FDR <= fdrThresh,
            , drop = F]
        enrich = enrichInfo
        if (!is.null(goTable)) {
            toDoI <- match(enrichInfo$term, goTable$term)
            enrichInfo$GOID = goTable$GOID[toDoI]
        }
        if (nrow(enrichInfo) == 0) {
            return()
        }
        enrichInfo$formattedLabel <- sapply(enrichInfo$term,
            function(curLabel) {
                curLabel <- strwrap(curLabel, labelWidth)
                paste(curLabel, collapse = "\n")
            })
        i <- match(enrichInfo$term, names(geneSets))
        if (any(is.na(i))) {
            stop("Could not find gene sets for ", sum(is.na(i)),
                " enriched terms.")
        }
        geneSets <- geneSets[i]
        if (is.null(enrichInfo$nGenes)) {
            enrichInfo$nGenes <- sapply(geneSets, length)
        }
        tmpSize <- -log10(enrichInfo$FDR)
        maxVal <- max(tmpSize[!is.infinite(tmpSize)])
        tmpSize[is.infinite(tmpSize)] <- maxVal + 2
        gsSizeRange <- range(tmpSize)
        if (gsSizeRange[1] == gsSizeRange[2]) {
            gsSizeRange[1] <- -log10(fdrThresh)
            gsSizeRange[2] <- gsSizeRange[2] + 1
        }
        tmpSize <- (tmpSize - gsSizeRange[1])/(gsSizeRange[2] -
            gsSizeRange[1])
        tmpSize <- nodeSizeRange[1] + tmpSize * (nodeSizeRange[2] -
            nodeSizeRange[1])
        enrichInfo$size <- round(tmpSize, 2)
        if (nrow(enrichInfo) == 1) {
            enrichInfo$cluster <- CLUST.COL[1]
            edgeMat <- NULL
        }
        else {
            pairI <- getUniquePairs(length(geneSets))
            distVal <- apply(pairI, 1, function(onePair) {
                overlapCoeff(geneSets[onePair])
            })
            distVal[distVal < overlapThresh] <- 0
            edgeMat <- data.frame(nodeA = pairI[, 1], nodeB = pairI[,
                2], coeff = distVal)
            enrichInfo$cluster <- prunedCol
            if (is.null(enrichInfo$pruneOutcome)) {
                termI <- 1:nrow(enrichInfo)
            }
            else {
                termI <- which(enrichInfo$pruneOutcome == enrichInfo$term)
            }
            if (length(termI) == 1) {
                enrichInfo$cluster[termI] <- CLUST.COL[1]
            }
            else {
                i <- which((edgeMat$nodeA %in% termI) & (edgeMat$nodeB %in%
                  termI))
                enrichInfo$id = termI
                g = igraph::graph_from_data_frame(edgeMat[which(edgeMat$coeff !=
                  0), ], directed = F, vertices = enrichInfo$id)
                adj = igraph::as_adjacency_matrix(g)
                clusters = igraph::clusters(g)
                clusters = split(names(clusters$membership),
                  clusters$membership)
                clusters <- lapply(clusters, as.numeric)
                lens <- sapply(clusters, length)
                clusters <- data.frame(id = unlist(clusters),
                  cluster = CLUST.COL[rep(1:length(clusters),
                    lens)], stringsAsFactors = F)
                enrichInfo$cluster[clusters$id] <- clusters$cluster
            }
            edgeMat <- edgeMat[edgeMat$coeff > 0, , drop = F]
            if (nrow(edgeMat) > 0) {
                edgeMat$size <- (edgeMat$coeff - overlapCoeffRange[1])/(overlapCoeffRange[2] -
                  overlapCoeffRange[1])
                edgeMat$size <- edgeWidthRange[1] + edgeMat$size *
                  (edgeWidthRange[2] - edgeWidthRange[1])
                edgeMat$coeff <- round(edgeMat$coeff, 2)
                edgeMat$size <- round(edgeMat$size, 2)
            }
            else {
                edgeMat <- NULL
            }
        }
        otherI <- order(enrichInfo$cluster)
        otherI <- otherI[order(enrichInfo$FDR[otherI])]
        termI <- which(enrichInfo$cluster[otherI] != prunedCol)
        if (length(termI) < length(otherI)) {
            otherI <- c(otherI[termI], otherI[-termI])
        }
        enrichInfo$id <- 1:nrow(enrichInfo)
        enrichInfo <- enrichInfo[otherI, , drop = F]
        enrichInfo$geneSetFraction <- round(enrichInfo$geneSetFraction *
            100, 1)
        enrichInfo$querySetFraction <- round(enrichInfo$querySetFraction *
            100, 1)
        if (!is.null(edgeMat)) {
            nam = c("source", "target", "label", "overlapCoeff",
                "width")
            orig = c("nodeA", "nodeB", "label", "coeff", "size")
            src = names(geneSets)[edgeMat$nodeA]
            trg = names(geneSets)[edgeMat$nodeB]
            edgeMat$label = paste(src, "(overlap)", trg)
            m = match(names(edgeMat), orig)
            names(edgeMat) = nam[m]
        }
        output = list(enrichInfo = enrichInfo, edgeMat = edgeMat)
        return(output)
    }

    fdrThresh = as.numeric(fdrThresh)
    bp_file = file.path(bp_path)
    queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index ==
        1)]))
    uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))
    if (!is.null(bp_input)) {
        bp = bp_input
    }
    else {
        bp <- readRDS(bp_file)
    }
    enrichMat.mn <- hyperG(querySet = queryGenes.mn, geneSets = bp,
        uni = uniGenes.mn, scoreMat = score, minSetSize = minSetSize,
        maxSetSize = maxSetSize, uniSize = NA)
    curr_exp = colnames(mat)[coln]
    queryGeneSets = list()
    queryGeneSets[[curr_exp]] = queryGenes.mn
    enrichMat.mn$filename <- colnames(mat)[coln]
    enrichMat_Ordered = enrichMat.mn[with(enrichMat.mn, order(FDR,
        -foldEnrichment)), ]
    scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),
        ]
    scoreMat <- scoreMat[match(uniGenes.mn, scoreMat$gene), "score",
        drop = F]
    rownames(scoreMat) <- uniGenes.mn
    colnames(scoreMat) <- curr_exp
    bp <- lapply(bp, intersect, uniGenes.mn)
    lens <- sapply(bp, length)
    bp <- bp[lens >= minSetSize & lens <= maxSetSize]
    go_file = file.path(go_path)
    if (!is.null(go_path))
        goTable = utils::read.delim(go_file, stringsAsFactors = F,
            check.names = F)
    if (!is.null(go_input)) {
        goTable = go_input
    }
    if (is.null(go_input) & is.null(go_path))
        goTable = NULL
    q = clusterEnrich(enrichInfo = enrichMat.mn, geneSets = bp,
        fdrThresh = fdrThresh, overlapThresh = 0.5, goTable = goTable)
    edgeMat = q$edgeMat
    enrichInfo = q$enrichInfo
    if (!is.null(enrichInfo)) {
        enrichInfo = enrichInfo %>% dplyr::arrange(FDR)
        enrichInfo$nOverlap = round(enrichInfo$geneSetFraction *
            enrichInfo$nGenes/100)
        enrichInfo$nQuery = round((enrichInfo$geneSetFraction *
            enrichInfo$nGenes/100)/(enrichInfo$querySetFraction/100))
        w = which(names(enrichInfo) %in% c("querySetFraction",
            "geneSetFraction", "foldEnrichment", "P", "FDR"))
        enrichInfo[, c("querySetFraction", "geneSetFraction",
            "foldEnrichment")] = round(enrichInfo[, c("querySetFraction",
            "geneSetFraction", "foldEnrichment")], 2)
        enrichInfo[, c("P", "FDR")] = signif(enrichInfo[, c("P",
            "FDR")], digits = 3)
    }
    else {
        print("no GO enrichment!")
    }
    nam = c("filename", "GOID", "term", "nGenes", "nQuery", "nOverlap",
        "querySetFraction", "geneSetFraction", "foldEnrichment",
        "P", "FDR", "overlapGenes", "maxOverlapGeneScore", "cluster",
        "id", "size", "formattedLabel")
    if (!is.null(goTable)) {
        m = match(enrichInfo$term, goTable$term)
        table(is.na(m))
        enrichInfo$GOID = goTable$GOID[m]
        enrichInfo = enrichInfo[, nam]
    }
    if (is.null(goTable)) {
        nam2 = c("filename", "term", "nGenes", "nQuery", "nOverlap",
            "querySetFraction", "geneSetFraction", "foldEnrichment",
            "P", "FDR", "overlapGenes", "maxOverlapGeneScore",
            "cluster", "id", "size", "formattedLabel")
        enrichInfo = enrichInfo[, nam2]
    }
    return(list(enrichInfo = enrichInfo, edgeMat = q$edgeMat))
}
