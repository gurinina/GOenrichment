# 'Computes Gene ontology enrichment using the hypergeometric test
#' The `hyperG` function evaluates a matrix of fitness scores to identify
#' genes with scores above a specified significance threshold.
#' @param querySet character vector of genes in query set, genes with significant fitness defect scores
#' @param geneSets A numeric index specifying the column of `mat` to analyze for significance.
#' @param scoreMat Dataframe of gene scores - first column = gene; column - second column = index, 1 and 0s designating significance - third column = fitness scores, descending order - can be NULL
#' see function `compSCORE`
#' @param maxSetSize Upper limit on the number of genes in a geneset included in the analysis (after restricting to the gene universe)
#' @param minSetSize Lower limit on the number of genes in a geneset included in the analysis (after restricting to the gene universe)
#' @param uni Character vector of genes in the universe (i.e. background set) - if NULL, must specify uniSize
#' @param uniSize The number of genes in the universe or length(uni)
#' @return Dataframe of enrichment results, sorted by increasing FDR value
#' if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
#' \describe{
#'   \item{querySetFraction}{the fraction of the query set that overlaps with the term set}
#'   \item{geneSetFraction}{the fraction of the term set that overlaps with the query set}
#'   \item{foldEnrichment}{the fold enrichment of the query set with the term genes}
#'   \item{P}{P value estimating the significance with which the query set is enriched with the term genes
#' FDR = FDR value estimating the significance of enrichment}
#'   \item{FDR}{FDR value estimating the significance of enrichment }
#'   \item{overlapGenes}{A |-separated list of genes in the overlap of the query set and the term set}
#' }
#' @importFrom dplyr arrange desc
#' @examples
#' # Example usage:
#' mat <- matrix(c(0.5, 1.2, 1.5, 0.8, 1.8, 0.3), nrow = 3, byrow = TRUE)
#' rownames(mat) <- c("Gene1", "Gene2", "Gene3")
#' scoreMat <- compSCORE(mat, coln = 2, sig = 1)
#' queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index >= 1)]))
#' uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))
#' enrichMat.mn <- hyperG(querySet = queryGenes.mn, geneSets = hGOBP.gmt,
#' @export
hyperG <- function (querySet, geneSets, uni, scoreMat,
            minSetSize = 5, maxSetSize = 300, uniSize = NA){
    if (!is.null(uni)) {
        geneSets <- lapply(geneSets, intersect, uni)
        lens <- sapply(geneSets, length)
        geneSets <- geneSets[lens >= minSetSize & lens <= maxSetSize]
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
        pVal <- phyper(length(overlapSet) - 1, length(geneSet),
            uniSize - length(geneSet), length(querySet), lower.tail = F)
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
    enrichCol <- data.frame(term = names(geneSets),
      querySetFraction = as.numeric(enrichInfo[,6]),
      geneSetFraction = as.numeric(enrichInfo[, 2]),
      foldEnrichment = as.numeric(enrichInfo[,3]),
      P = as.numeric(enrichInfo[, 4]),
      FDR = p.adjust(as.numeric(enrichInfo[, 4]), method = "BH"),
      overlapGenes = enrichInfo[, 1],
      maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = F)
    rownames(enrichCol) <- NULL
    enrichCol = enrichCol[order(enrichCol$FDR), ]
}
