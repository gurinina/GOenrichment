



#' @title runGORESP
#' @description This function performs a GO enrichment analysis using the hypergeometric
#' test for the set of query genes passing the user-defined fitness score
#' threshold compared to the genes in the background set (universe).
#' @param mat numeric matrix of fitness data
#' @param coln mumeric. the column of the matrix with sample of interest
#' @param curr_exp character, name of exp, usually after the sample name so set to colnames(mat)[coln]
#' @param sig mumeric, significance threshold of fitness defect score, default is 1
#' @param fdrThresh numeric, fdr threshold for significance cutoff of enrichments, default = 0.2
#' @param bp_path character, path of gmt file
#' @param bp_input gmt file or NULL, **Note: bp_path and bp_input can't both be NULL
#' @param go_path character, path of file with GOID and terms path of gmt file alternatively bp_input, gmt file itself,
#' @param go_input dataframe with GOID and terms matching bp_input  **Note: go_path and go_input can both be NULL
#' @param minGeneSetSize numeric, lower limit on the number of genes in a geneset included in the analysis, default = 5
#' @param maxGeneSetSize numeric, upper limit on the number of genes in a geneset included in the analysis,default = 300
#' @return RETURNS two dataframes: enrichInfo with enrichment results and edgeMat with information about
#' geneset term overlap and clusters.
#' dataframe includes the columns: filename,term, nGenes, nQuery, nOverlap,
#'v@details columns enrichInfo:
#'     filename
#'     GOID: GO ID
#'     term: GO term
#'     nGenes: # genes in the geneset
#'     nQuery: # genes in the query
#'     nOverlap: # genes overlapping the query and the geneset
#'     querySetFraction: the fraction of the query set (the set that passes the significance threshold)
#'      that overlaps with the term set
#'     geneSetFraction: the fraction of the term set that overlaps with the query set
#'     foldEnrichment = the fold enrichment of the query set over the bgRate, where bgRate
#'        is nGenes in geneSet/genes in the universe
#'     P = P value estimating the significance with which the query set is enriched with the term genes
#'     FDR = FDR value estimating the significance of enrichment
#'     overlapGenes = a |-separated list of genes that overlap the query set and the geneSet
#'                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#'	   maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
#'	   cluster: cluster color based on overlap of term ID, set as 50%
#'	   id: cluster ID
#'	   size: node size for downstream plotting based on FDR score
#'	   formatted label: for downstream plotting
#' columns edgeMat:
#'	   source
#'	   target
#'	   overlapCoeff
#'	   width
#'	   label
#'	   These map GO terms from the enrichInfo serving to show the
#'	   overlap between terms for downstream visualization purposes
#'	   and for understanding the relationship between terms
#'v@export
runGORESP = function (mat,coln,curr_exp = colnames(mat)[coln], sig = 1, fdrThresh = 0.2, bp_path = NULL,bp_input = NULL,
                      go_path = NULL,go_input = NULL,minGeneSetSize = 5,maxGeneSetSize = 300){


##############
CLUST.COL <- c("#FF00CC","#33CCFF", "#33CC00", "#9900FF", "#FF9900", "#FFFF00", "#FFCCFF", "#FF0000", "#006600", "#009999", "#CCCC00", "#993300", "#CC99CC", "#6699CC","#CCCCFF", "#FFCC99", "#9966FF", "#CC6600", "#CCFFFF", "#99CC00", "#FF99FF", "#0066FF", "#66FFCC", "#99CCFF", "#9999CC", "#CC9900", "#CC33FF", "#006699", "#F5DF16", "#B5185E", "#99FF00", "#00FFFF", "#990000", "#CC0000", "#33CCCC", "#CC6666", "#996600", "#9999FF", "#3366FF")

prunedCol <- "#BEBEBE"





##### required functions
# computes the number of unique pairs given the number of items to consider
# maxVal - the maximum number of items
# RETURNS a 2-column matrix where each row contains a different pair, specified with item indices
getUniquePairs = function (maxVal)
{
  firstI <- rep(1:(maxVal - 1), (maxVal - 1):1)
  secondI <- sapply(2:maxVal, function(x) {
    x:maxVal
  })
  cbind(firstI, unlist(secondI))
}



##############


# computes enrichment using the hypergeometric test, and uses the resulting P values with
# the Benjamini Hochberg method to estimate FDR values
# querySet - character vector of genes in query set
# geneSets - named list of gene sets to test for significant overlap w/ the query set
# scoreMat - dataframe of gene scores
#          - first column = scores, gene column
#          - can be NULL
# uni - character vector of genes in the universe (i.e. background set)
#     - if NULL, must specify uniSize
# uniSize - the # of genes in the universe
# minSetSize, maxSetSize - min/max # of genes in geneSets (after restricting to the gene universe)
# RETURNS a dataframe of enrichment results, sorted by increasing FDR value. The columns are:
#         term = name of gene set
#         querySetFraction = the fraction of the query set that overlaps with the term set
#         geneSetFraction = the fraction of the term set that overlaps with the query set
#         foldEnrichment = the fold enrichment of the query set with the term genes
#         P = P value estimating the significance with which the query set is enriched with the term genes
#         FDR = FDR value estimating the significance of enrichment
#         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#	  maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
hyperG = function (querySet, geneSets, uni, scoreMat, minSetSize = minGeneSetSize,
  maxSetSize = maxGeneSetSize, uniSize = NA)
{
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
    pVal <- stats::phyper(length(overlapSet) - 1, length(geneSet),
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
  enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[,
    6]), geneSetFraction = as.numeric(enrichInfo[, 2]), foldEnrichment = as.numeric(enrichInfo[,
      3]), P = as.numeric(enrichInfo[, 4]), FDR = stats::p.adjust(as.numeric(enrichInfo[,
        4]), method = "BH"), overlapGenes = enrichInfo[, 1],
    maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = F)
  rownames(enrichCol) <- NULL
  enrichCol = enrichCol[order(enrichCol$FDR), ]
}
#



#######the overlap of genesets for all combinations
#######If set X is a subset of Y or the converse then the overlap coefficient is equal to 1.

# compute the overlap coefficient given a pair of (gene) sets
# gsPairList - a list of two sets (each set is a vector of IDs)
# RETURNS the overlap coefficient

overlapCoeff = function (gsPairList)
{
  length(intersect(gsPairList[[1]], gsPairList[[2]]))/min(length(gsPairList[[1]]),
    length(gsPairList[[2]]))
}

####### given enrichInfo after hyperG, computes edgeMat for making enrichment map
# generates an enrichment map in xgmml format using hypergeometric test statistics
# enrichInfo - dataframe with enrichment stats for gene sets (one per row), with the following columns:
#              term, geneSetFraction, querySetFraction, FDR, overlapGenes, maxOverlapGeneScore
#            - see documentation for the output of hyperG() for descriptions of these columns
# geneSets - named list of gene sets tested for significant overlap w/ the query set,
#              restricted to genes in the universe
#
# fdrThresh - FDR threshold; only show gene sets that pass this significance threshold
# overlapThresh - an edge between a pair of enriched gene sets will only be shown if the overlap coefficient
#              is >= overlapThresh



# goTable - a dataframe with the following columns describing GO terms:
#         - "term" (GO term), "id" (GOID)
#         - if provided (i.e. not NULL), the GO ID numbers of the enriched GO terms will be saved in
#           the output xgmml file as a node attribute called "GOID".
#         - the GOID allows for an easy link to the GO website page for the associated GO term

###############
#### query set is genes with significant fitness defect, uni is all the genes in the data matrix, the union
clusterEnrich = function (enrichInfo, geneSets, fdrThresh = 0.1, overlapThresh = 0.5,
                            go_path = go_path, go_input = go_input){


  go_file = file.path(go_path)
  if(!is.null(go_path)) goTable = utils::read.delim(go_file,stringsAsFactors = F,check.names = F)
  if(!is.null(go_input))  {goTable = go_input}


  nodeSizeRange <- c(10, 40)
  prunedCol <- "#BEBEBE"
  labelWidth <- 20
  edgeWidthRange <- c(1, 5)
  overlapCoeffRange <- c(overlapThresh, 1)



  enrichInfo <- enrichInfo[enrichInfo$FDR <= fdrThresh, , drop = F]

  enrich = enrichInfo


  if(!is.null(goTable)) {
    toDoI <- match(enrichInfo$term, goTable$term)
    enrichInfo$GOID = goTable$GOID[toDoI]
  }

  if (nrow(enrichInfo) == 0) {

    return()
  }
  enrichInfo$formattedLabel <- sapply(enrichInfo$term, function(curLabel) {
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
  tmpSize <- (tmpSize - gsSizeRange[1])/(gsSizeRange[2] - gsSizeRange[1])
  tmpSize <- nodeSizeRange[1] + tmpSize * (nodeSizeRange[2] -nodeSizeRange[1])
  enrichInfo$size <- round(tmpSize, 2)
  if (nrow(enrichInfo) == 1) {
    enrichInfo$cluster <- CLUST.COL[1]
    edgeMat <- NULL
  }



#########################   EDGE MAT   ########################
  else {
    pairI <- getUniquePairs(length(geneSets))
    distVal <- apply(pairI, 1, function(onePair) {
      overlapCoeff(geneSets[onePair])
    })
    distVal[distVal < overlapThresh] <- 0
    edgeMat <- data.frame(nodeA = pairI[, 1], nodeB = pairI[,2], coeff = distVal)
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

      g=igraph::graph_from_data_frame(edgeMat[which(edgeMat$coeff!=0),],directed = F,vertices = enrichInfo$id)
      adj = igraph::as_adjacency_matrix(g)
      clusters = igraph::clusters(g)
      clusters = split(names(clusters$membership),clusters$membership)


      clusters <- lapply(clusters, as.numeric)

      lens <- sapply(clusters, length)
      clusters <- data.frame(id = unlist(clusters), cluster = CLUST.COL[rep(1:length(clusters),
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


  #print("edgeMat is NULL")


  if (!is.null(edgeMat)) {
    nam = c("source","target","label","overlapCoeff","width")
    orig = c("nodeA","nodeB","label","coeff","size")

    src = names(geneSets)[edgeMat$nodeA]

    trg = names(geneSets)[edgeMat$nodeB]
    edgeMat$label = paste(src,"(overlap)",trg)
    m = match(names(edgeMat),orig)
    names(edgeMat) = nam[m]
  }


  output = list(enrichInfo = enrichInfo,edgeMat = edgeMat)

  return(output)
}
###############

compSCORE <- function(mat,coln, sig = 1){

  df = data.frame(score = mat[,coln],stringsAsFactors = F)

  df$gene = rownames(mat)
  rownames(df) = df$gene
  df$index=0
  wdf = which(df$score >= sig)
  df$index[wdf]=1
  df = df[,c('index','score','gene')]
  df = dplyr::arrange(df,dplyr::desc(score))
  df
}

  ###### bulk of code

  score = compSCORE(mat,coln,sig = sig)

  fdrThresh = as.numeric(fdrThresh)

  bp_file = file.path(bp_path)
  scoreMat = score


  queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index == 1)]))
  uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))

  if(!is.null(bp_input))  {bp = bp_input} else {bp <- readRDS(bp_file)}

  #uniGenes.mn <- unique(intersect(scoreMat$gene,unlist(bp,use.names=F)))

  #### intersect of the geneSets with the backgroundSet, filtering for size


  enrichMat.mn <- hyperG(querySet = queryGenes.mn, geneSets = bp,
                           uni = uniGenes.mn, scoreMat = score, minSetSize = minGeneSetSize,
                           maxSetSize = maxGeneSetSize, uniSize = NA)
  curr_exp = colnames(mat)[coln]
  queryGeneSets = list()
  queryGeneSets[[curr_exp]] = queryGenes.mn
  enrichMat.mn$filename <- colnames(mat)[coln]
  enrichMat_Ordered = enrichMat.mn[with(enrichMat.mn, order(FDR,
                                                            -foldEnrichment)), ]
  scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),]
  scoreMat <- scoreMat[match(uniGenes.mn, scoreMat$gene), "score",drop = F]
  rownames(scoreMat) <-   uniGenes.mn
  colnames(scoreMat) <- curr_exp
#
#   nonEnrichMat.mn <- genesNotInEnrichedTerm(queryGeneSets,
#                                               enrichMat.mn, scoreMat, NONSPECIFIC.TERMS$bp, fdrThresh)
#
  maxSetSize = maxGeneSetSize
  #### intersect of the geneSets with the backgroundSet, filtering for size
  bp <- lapply(bp, intersect, uniGenes.mn)
  lens <- sapply(bp, length)
  bp <- bp[lens >= minGeneSetSize & lens <= maxGeneSetSize]

  q = clusterEnrich(enrichInfo = enrichMat.mn, geneSets = bp,
                      fdrThresh = fdrThresh, overlapThresh = 0.5,
                      go_path = go_path,go_input = go_input)

  edgeMat = q$edgeMat
  enrichInfo = q$enrichInfo
  m = match(enrichInfo$term,go_input$term)
      table(is.na(m))
      enrichInfo$GOID = go_input$GOID[m]

  if(!is.null(enrichInfo)) {

    enrichInfo = enrichInfo %>% dplyr::arrange(FDR)
    enrichInfo$nOverlap = round(enrichInfo$geneSetFraction*enrichInfo$nGenes/100)
    enrichInfo$nQuery = round((enrichInfo$geneSetFraction*enrichInfo$nGenes/100)/(enrichInfo$querySetFraction/100))

    w = which(names(enrichInfo) %in% c("querySetFraction", "geneSetFraction" ,
                                   "foldEnrichment" , "P" , "FDR" ))
    enrichInfo[,c("querySetFraction","geneSetFraction", "foldEnrichment")] =
      round(enrichInfo[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2)

    enrichInfo[,c("P", "FDR")] =
      signif(enrichInfo[,c("P", "FDR")],digits = 3)

  }else { print("no GO enrichment!") }
    nam = c(
      "filename",
      "GOID",
      "term",
      "nGenes"  ,
      "nQuery",
      "nOverlap"     ,
      "querySetFraction",
      "geneSetFraction" ,
      "foldEnrichment"  ,
      "P"            ,
      "FDR"      ,
      "overlapGenes",
      "maxOverlapGeneScore" ,
      "cluster"      ,
      "id" ,
      "size"    ,
      "formattedLabel"
    )

    go_file = file.path(go_path)
    if(!is.null(go_path)) goTable = utils::read.delim(go_file,stringsAsFactors = F,check.names = F)
    if(!is.null(go_input))  {goTable = go_input}

    # wnam = which(names(enrichInfo)%in% nam)
    # print(nam[wnam])


    if(!is.null(goTable)){

      m = match(enrichInfo$term,goTable$term)
      table(is.na(m))
      enrichInfo$GOID = goTable$GOID[m]


    enrichInfo = enrichInfo[,nam]}
    sdf = setdiff(nam,names(enrichInfo))

    if(is.null(go_path)&is.null(go_input)) {nam2 = c(
      "filename",

      "term",
      "nGenes"  ,
      "nQuery",
      "nOverlap"     ,
      "querySetFraction",
      "geneSetFraction" ,
      "foldEnrichment"  ,
      "P"            ,
      "FDR"      ,
      "overlapGenes",
      "maxOverlapGeneScore" ,
      "cluster"      ,
      "id" ,
      "size"    ,
      "formattedLabel"
    )

     enrichInfo = enrichInfo[,nam2]}

  return(list(enrichInfo = enrichInfo , edgeMat = q$edgeMat))
}


#' @title hyperG
#' @description computes enrichment using the hypergeometric test given:
#' @param querySet character vector of genes in query set, genes with significant fitness defect scores
#' @param geneSets named list of gene sets to test for significant overlap w/ the query set
#' @param scoreMat dataframe of gene scores
#'          - first column = gene column
#'          - second column = index, 1 and 0s designating significance
#'          - third column = fitness scores, descending order
#'          - can be NULL
#' @param uni character vector of genes in the universe (i.e. background set)
#'     - if NULL, must specify uniSize
#' @param uniSize the number of genes in the universe
#' @param maxSetSize upper limit on the number of genes in a geneset included in the analysis
#' (after restricting to the gene universe)
#' @param minSetSize lower limit on the number of genes in a geneset included in the analysis
#' (after restricting to the gene universe)
#'
#' @return dataframe of enrichment results, sorted by increasing FDR value.
#' @details The columns of the dataframe are:
#          term = name of gene set
#'         querySetFraction = the fraction of the query set that overlaps with the term set
#'         geneSetFraction = the fraction of the term set that overlaps with the query set
#'         foldEnrichment = the fold enrichment of the query set with the term genes
#'         P = P value estimating the significance with which the query set is enriched with the term genes
#'         FDR = FDR value estimating the significance of enrichment
#'         overlapGenes = a |-separated list of genes in the overlap of the query set and the term set;
#'                        if scoreMat is provided (not NULL), the scores of the genes are shown in parentheses
#'	       maxOverlapGeneScore = if scoreMat is provided (not NULL), the maximum score of the overlapGenes
#' @export
hyperG = function (querySet, geneSets, uni, scoreMat, minSetSize = 5,
  maxSetSize = 300, uniSize = NA)
{
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
  enrichCol <- data.frame(term = names(geneSets), querySetFraction = as.numeric(enrichInfo[,
    6]), geneSetFraction = as.numeric(enrichInfo[, 2]), foldEnrichment = as.numeric(enrichInfo[,
      3]), P = as.numeric(enrichInfo[, 4]), FDR = p.adjust(as.numeric(enrichInfo[,
        4]), method = "BH"), overlapGenes = enrichInfo[, 1],
    maxOverlapGeneScore = as.numeric(enrichInfo[, 5]), stringsAsFactors = F)
  rownames(enrichCol) <- NULL
  enrichCol = enrichCol[order(enrichCol$FDR), ]
}

#'@title compSCORE:
#'@description computes a dataframe of designating significant fitness scores from a
#'matrix of screening data
#'@param mat numeric matrix of fitness scores
#'@param coln numeric, column # of matrix with the sample of interest
#'@param sig numeric, significance threshold, default = 1
#'@return dataframe with 3 columns: gene, index, score. index indicating significant 1,
#' or nonsignificant, 0
#'@export
#'@importFrom dplyr "%>%"

compSCORE <- function(mat,coln, sig = 1){
   df = data.frame(score = mat[,coln],stringsAsFactors = F)

   df$gene = rownames(mat)
   rownames(df) = df$gene
   df$index=0
   wdf = which(df$score >= sig)
   df$index[wdf]=1
   df = df[,c('index','score','gene')]
   df = df %>% dplyr::arrange(dplyr::desc(score))
   df
   }
############
#' @title visSetup
#' @description this function takes the output from runGORESP generates objects for drawing GO enrichment networks
#' @param enrichInfo a data frame, see ??runGORESP enrichInfo output for explanation of these columns
#' @param edgeMat a data frame, see ??runGORESP edgeMat output for explanation of these columns
#' @param fontsize fontsize
#' @param fontface fontface
#' @return returns two dataframes in a list, vis$nodes with the following columns:
#' id, filename,term, nGenes, nQuery, nOverlap, querySetFraction, geneSetFraction, foldEnrichment,
#' P, FDR, overlapGenes, maxOverlapGeneScore: see ??runGORESP enrichInfo output for explanation of these columns
#' visualization columns: size, formattedLabel, label, color.background,color.border,
#' borderWidthSelected, color.highlight.border, color.highlight.background,
#' color.hover.background, color.hover.border, font.face, shape, font.size, font.bold,
#' borderWidth, labelHighlightBold and and vis$edges with the columns from, to, overlapCoeff,
#' width and color
#' @export
visSetup = function(enrichInfo, edgeMat, fontsize = 22, fontface = "Arial") {


  n = enrichInfo
  e = edgeMat

  w = which(names(n) == "id")
  coln = (1:ncol(n))[-w]
  n = n[, c(w, coln)]

  if (is.null(e) & !is.null(n)) {
    gr = igraph::make_empty_graph(nrow(enrichInfo))
    v = gr
    igraph::V(v)$color.background = n$cluster
    v = igraph::set_vertex_attr(v, "label", value = n$formattedLabel)

  }


  if (!is.null(e) & !is.null(n)) {
    w = which(names(e) == "label")
    let = igraph::graph_from_data_frame(e[, -w], vertices = n, directed = F)
    v = igraph::set_vertex_attr(let, "label", value = n$formattedLabel)
    igraph::V(v)$color.background = n$cluster
  }

  vis = visNetwork::toVisNetworkData(v)
  vis$nodes = data.frame(vis$nodes, stringsAsFactors = F)
  if (!is.null(vis$edges))
    vis$edges = data.frame(vis$edges, stringsAsFactors = F)
  m = match(vis$nodes$id, n$id)
  vis$nodes$label = n$formattedLabel[m]
  vis$nodes$FDR = n$FDR[m]
  vis$nodes$FDR = signif(vis$nodes$FDR, 5)


  w = which(duplicated(vis$nodes$FDR))
  if (length(w) > 0) {
    vis$nodes =  dplyr::group_by(vis$nodes,FDR)
    vis$nodes = dplyr::mutate(vis$nodes, jitter = if (n() > 1) abs(jitter(FDR)) else FDR)
    w = which(names(vis$nodes) == "FDR")
    vis$nodes = vis$nodes[, -w]
    w = which(names(vis$nodes) == "jitter")
    names(vis$nodes)[w] = "FDR"
    vis$nodes$FDR = signif(vis$nodes$FDR, 6)
  }

  w = which(duplicated(vis$nodes$FDR))



  vis$nodes$term = n$term[m]
  vis$nodes$nOverlap = n$nOverlap[m]
  vis$nodes$nQuery= n$nQuery[m]
  vis$nodes$nGenes = n$nGenes[m]
  vis$nodes$geneSetFraction = n$geneSetFraction[m]
  vis$nodes$querySetFraction = n$querySetFraction[m]
  vis$nodes$filename =  n$filename[m]

  vis$nodes$formattedLabel = n$formattedLabel[m]
  vis$nodes$overlapGenes = n$overlapGenes[m]
  vis$nodes$label = vis$nodes$formattedLabel
  vis$nodes$color.border = "black"
  vis$nodes = dplyr::arrange(vis$nodes,label)
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
  vis$nodes = dplyr::arrange(vis$nodes,FDR)
  vis

}
############

#' @title runNetwork
#' @description Uses the output from visSetup to generate a network of gene enrrichments.
#' @param nodes dataframe returned from visSetup, "nodes" see ?? visSetup
#' @param edges dataframe returned from visSetup, "edges" see ?? visSetup
#' @param height A numeric value
#' @return network
#' @details The larger a node in the netwark, the more significance the
#' GO enrichment for that term. The thicker the edge connecting twa nodes,
#' the greater the overlap betweem twp terms. Bold and enlarged font size
#' indicate the most enriched term for a given cluster.
#' @export
runNetwork <- function(nodes,edges, height = 1300, main = list(text = nodes$filename[1],
style = "font-size:40px;text-align:center;"),...){
n = nodes

n <- n%>% dplyr::arrange(term)
names=n$id

visNetwork(nodes, edges, width = "100%",height = height, main = main,...) %>%
visNodes(shadow=list(enabled=T,size=25)) %>%

visOptions(
highlightNearest = list(enabled = T, degree = 5, hover = T),
nodesIdSelection = list(enabled = TRUE, values = names,
style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),
selectedBy = list(variable="FDR",
style = 'width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%
visIgraphLayout(type = "full")  %>%
visEvents(select = "function(nodes) {
Shiny.onInputChange('current_node_id', nodes.nodes);
;}")
}




