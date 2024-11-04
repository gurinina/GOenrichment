#' Computes Gene Ontology enrichment using the hypergeometric test
#'
#' The `hyperG` function identifies enriched gene sets by comparing a query set of genes with significant fitness defect scores to a background universe of genes. The hypergeometric test assesses the significance of overlap between the query set and predefined gene sets.
#'
#' @param querySet Character vector of genes in the query set with significant fitness scores.
#' @param geneSets List of gene sets (pathways or processes), each represented as a character vector of gene symbols.
#' @param scoreMat Data frame of gene scores. The first column is `gene`, the second column (`index`) indicates significance (1 for significant, 0 otherwise), and the third column is `score`, representing fitness scores in descending order. Can be NULL.
#' @param maxSetSize Numeric, upper limit on the number of genes in a gene set included in the analysis after restricting to the universe.
#' @param minSetSize Numeric, lower limit on the number of genes in a gene set included in the analysis after restricting to the universe.
#' @param uni Character vector of genes in the universe (background set). If NULL, `uniSize` must be specified.
#' @param uniSize Numeric, the number of genes in the universe or length of `uni` if `uni` is provided.
#' @return Data frame of enrichment results, sorted by increasing FDR value.
#' \describe{
#'   \item{term}{The name of the gene set.}
#'   \item{querySetFraction}{Fraction of genes in the query set that are also in the gene set.}
#'   \item{geneSetFraction}{Fraction of the gene set that overlaps with the query set.}
#'   \item{foldEnrichment}{Fold enrichment of the query set with the term genes.}
#'   \item{P}{P-value estimating the significance with which the query set is enriched with the term genes.}
#'   \item{FDR}{False discovery rate (adjusted P-value) for the enrichment significance.}
#'   \item{overlapGenes}{`|-` separated list of genes in the overlap between the query set and the term set.}
#'   \item{maxOverlapGeneScore}{Maximum fitness score among the overlapping genes.}
#' }
#' @importFrom dplyr arrange desc
#' @examples
#' queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index == 1)]))
#' uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))
#' enrichMat.mn <- hyperG(querySet = queryGenes.mn, geneSets = hGOBP.gmt, uni = uniGenes.mn, scoreMat = score, minSetSize = minSetSize)
#' @export
hyperG <- function(querySet, geneSets, uni, scoreMat, minSetSize = 5, maxSetSize = 300, uniSize = NA) {
    if (!is.null(uni)) {
        geneSets <- lapply(geneSets, intersect, uni)
        lens <- sapply(geneSets, length)
        geneSets <- geneSets[lens >= minSetSize & lens <= maxSetSize]
        uniSize <- length(uni)
    }
    if (!is.null(scoreMat)) {
        scoreMat <- scoreMat[order(scoreMat$score, decreasing = TRUE), ]
        if (!is.null(uni)) {
            i <- match(uni, scoreMat$gene)
            scoreMat <- scoreMat[sort(i[!is.na(i)]), ]
        }
        scoreMat$score <- round(scoreMat$score, 2)
    }
    enrichInfo <- sapply(geneSets, function(geneSet) {
        overlapSet <- intersect(querySet, geneSet)
        pVal <- phyper(length(overlapSet) - 1, length(geneSet), uniSize - length(geneSet), length(querySet), lower.tail = FALSE)
        if (length(overlapSet) > 0) {
            overlapSet <- sort(overlapSet)
        }
        overlapSize <- length(overlapSet)
        maxScore <- if (is.null(scoreMat)) NA else scoreMat$score[sort(match(overlapSet, scoreMat$gene))][1]
        overlapSet <- if (!is.null(scoreMat)) paste(paste(scoreMat$gene[match(overlapSet, scoreMat$gene)], "(", scoreMat$score[match(overlapSet, scoreMat$gene)], ")", sep = ""), collapse = "|") else paste(overlapSet, collapse = "|")
        bgRate <- length(geneSet) / uniSize
        foldEnrich <- overlapSize / length(querySet) / bgRate
        c(overlapSet, overlapSize / length(geneSet), foldEnrich, pVal, maxScore, overlapSize / length(querySet))
    })
    enrichInfo <- t(enrichInfo)
    enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[, 6]), geneSetFraction = as.numeric(enrichInfo[, 2]), foldEnrichment = as.numeric(enrichInfo[, 3]), P = as.numeric(enrichInfo[, 4]), FDR = p.adjust(as.numeric(enrichInfo[, 4]), method = "BH"), overlapGenes = enrichInfo[, 1], maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = FALSE)
    rownames(enrichCol) <- NULL
    enrichCol <- enrichCol[order(enrichCol$FDR), ]
}
