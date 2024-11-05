#' Calculate Significant Fitness Scores
#' @import dplyr
#' The `compSCORE` function evaluates a matrix of fitness scores to identify
#' genes with scores above a specified significance threshold.
#'
#' @param mat A numeric matrix where rows represent genes and columns represent different conditions or samples.
#' @param coln A numeric index specifying the column of `mat` to analyze for significance.
#' @param sig A numeric significance threshold; scores equal to or above this threshold are marked as significant (default is 1).
#' @return A data frame with three columns:
#' \describe{
#'   \item{index}{A binary indicator of significance: 1 if the score is significant, 0 otherwise.}
#'   \item{score}{The fitness score of the gene in the specified column.}
#'   \item{gene}{The name of the gene, taken from the row names of `mat`.}
#' }
#' @examples
#' # Example usage:
#' mat <- matrix(c(0.5, 1.2, 1.5, 0.8, 1.8, 0.3), nrow = 3, byrow = TRUE)
#' rownames(mat) <- c("Gene1", "Gene2", "Gene3")
#' compSCORE(mat, coln = 2, sig = 1)
#' @export
compSCORE <- function(mat, coln, sig = 1) {
    # Create a data frame with fitness scores
    df <- data.frame(score = mat[, coln], stringsAsFactors = FALSE)

    # Add gene names from rownames of mat
    df$gene <- rownames(mat)
    rownames(df) <- df$gene

    # Initialize significance index to 0 (not significant)
    df$index <- 0

    # Identify significant scores
    wdf <- which(df$score >= sig)
    df$index[wdf] <- 1

    # Arrange and return the data frame with index, score, and gene
    df <- df[, c("index", "score", "gene")]
    df <- df %>% dplyr::arrange(dplyr::desc(score))
    df
}
