#' Count number of times a term is annotated within gene sample.
#' @description To count the number of times a term occurs in a given sample
#' @description  Precursor of hypergeometric distribution calculation.
#' @param genes vector of hgnc gene names
#' @param annots output of annotateGO()
#'
#' @importFrom dplyr select
#' @export
#' @return dataframe of GO term counts
#'
GOdistr <- function(genes, annots = annotations){

  distr <- annots %>% dplyr::select(SYMBOL, TERM) %>%
    filter(SYMBOL %in% genes) %>%
    group_by(TERM) %>%
    summarise(SYMBOLS = paste(SYMBOL, collapse=", ")) %>%
    arrange(TERM)

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
