#' Homo sapiens GO Biological Process (BP) Gene Sets
#'
#' BP gene sets for Homo sapiens downloaded from Gary Bader's lab in July 2022.
#' This dataset excludes terms annotated by electronic annotation (IEA).
#'
#' @name hGOBP.gmt
#' @docType data
#' @title Homo sapiens GO Biological Process (BP) Genesets
#' @description This data provides GO biological process (BP) gene sets, where
#' each set corresponds to a BP term with genes annotated to it, minus those
#' with IEA annotations.
#' @format A list of BP gene sets where each name represents a GO term.
#' @usage data(hGOBP.gmt)
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
#' @examples
#' data(hGOBP.gmt)
#' head(hGOBP.gmt)
