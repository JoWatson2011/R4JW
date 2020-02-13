#' Count number of times a term is annotated within gene sample.
#' @description To count the number of times a term occurs in a given sample
#' @description  Precursor of hypergeometric distribution calculation.
#' @param genes vector of hgnc gene names
#' @param annots output of annotateGO()
#'
#' @importFrom dplyr select group_by summarise arrange
#' @importFrom rlang .data
#' @export
#' @return dataframe of GO term counts
#'
GOdistr <- function(genes, annots){
  distr <- annots %>% dplyr::select(.data$SYMBOL, .data$TERM) %>%
    filter(.data$SYMBOL %in% genes) %>%
    group_by(.data$TERM) %>%
    summarise(SYMBOLS = paste(.data$SYMBOL, collapse=", ")) %>%
    arrange(.data$TERM)

  distr$SYMBOLS <- apply(distr, 1, function(x){
    t <- strsplit(x["SYMBOLS"], ", ")
    x["SYMBOLS"] <- paste(unique(t[[1]]), collapse = ", ")
  })

  distr$COUNT <- apply(distr, 1, function(x){
    t <- strsplit(x["SYMBOLS"], ", ")
    x["COUNT"] <- length(t[[1]])
  })

  return(distr)
}
