#' GO enrichment analysis using the hypergeometric test
#' @importFrom magrittr %>%
#' Performs a GO enrichment analysis using the hypergeometric test for a set of query genes that pass a user-defined fitness score threshold. It compares the query genes against a background set (universe) and outputs enrichment and clustering information.
#'
#' @param scoreMat Data frame with gene scores. The first column is "gene," the second is "index" (indicating significance with 1 or 0), and the third is "score" (fitness scores in descending order).
#' @param curr_exp Character; experiment label for the analysis, default = "test".library(GOenrichment)
#' @param fdrThresh Numeric; FDR threshold for enrichment significance, default = 0.2.
#' @param bp_path Character; path to a .gmt file for biological process (BP) gene sets.
#' @param bp_input List; gene sets in .gmt format. If provided, `bp_path` is ignored.
#' @param go_path Character; path to file with GO terms and IDs, or NULL if using `go_input`.
#' @param go_input Data frame; GO terms and IDs. If provided, `go_path` is ignored.
#' @param minSetSize Numeric; minimum number of genes in a gene set for inclusion, default = 5.
#' @param maxSetSize Numeric; maximum number of genes in a gene set for inclusion, default = 300.
#' @param uniSize Numeric; number of genes in the universe, or length of `uni`.
#' @return List containing:
#' \describe{
#'   \item{enrichInfo}{Data frame with GO enrichment results including columns:}
#'   \describe{
#'     \item{term}{Name of the GO term.}
#'     \item{querySetFraction}{Fraction of the query set overlapping with the term set.}
#'     \item{geneSetFraction}{Fraction of the term set overlapping with the query set.}
#'     \item{foldEnrichment}{Fold enrichment of query genes within the GO term genes.}
#'     \item{P}{P-value for enrichment significance.}
#'     \item{FDR}{False discovery rate (adjusted P-value) for enrichment significance.}
#'     \item{overlapGenes}{List of overlapping genes in the query and term sets.}
#'     \item{maxOverlapGeneScore}{Maximum fitness score among overlapping genes, if `scoreMat` is provided.}
#'     \item{cluster, id, size, formattedLabel}{Additional clustering and visualization information.}
#'   }
#'   \item{edgeMat}{Data frame for gene set term overlaps, including source, target, and overlap coefficients for visualization.}
#' }
#' @examples
#' # Example usage:
#' scoreMat <- compSCORE(mat, coln, sig = sig)
#' curr_exp <- colnames(mat)[coln]
#' hresp <- GOenrichment::runGORESP(fdrThresh = 0.2, scoreMat = scoreMat,
#'   coln = 1, curr_exp = "test", sig = 1,
#'   bp_input = hGOBP.gmt, go_input = NULL, minSetSize = 20,
#'   maxSetSize = 300)
#' @export
runGORESP <- function (scoreMat, curr_exp = "test",
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
        uni = uniGenes.mn, scoreMat = scoreMat, minSetSize = minSetSize,
        maxSetSize = maxSetSize, uniSize = NA)
    curr_exp = "test"
    queryGeneSets = list()
    queryGeneSets[[curr_exp]] = queryGenes.mn
    enrichMat.mn$filename <- curr_exp
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
