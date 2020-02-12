#' Over Representation Analysis
#' @description Calculate over enrichment of GO terms in
#' @description a sample Explanation can be found in:
#' @description https://mengnote.blogspot.co.uk/2012/12/calculate-correct-hypergeometric-p.html
#' @param GO.clust
#' @param sample
#' @param HeLa.bkgd
#' @param HeLa.GO
#' @param p
#'
#' @return
#' @export
GOenrichment <- function(GO.clust, sample, HeLa.bkgd = 12427, HeLa.GO = GO.distr.HeLa, p = 0.05){

  suppressPackageStartupMessages(library("dplyr"))

  successes.sample <- GO.clust$COUNT    # How many times term 'i' occurred in sample

  successes.bkgd <- HeLa.GO[HeLa.GO$TERM %in% GO.clust$TERM,] %>% arrange(TERM)
  successes.bkgd <- successes.bkgd$COUNT   #How many times term 'i' occured in background

  failure <- HeLa.bkgd - successes.bkgd      #Number of peptides in background minus the peptides with that term

  sample.size <- nrow(sample)   # Number of peptides in sample

  GO.clust$p <- phyper(successes.sample-1,
                       successes.bkgd,
                       failure,
                       sample.size,
                       lower.tail = F)

  GO.clust$adjustedP <- p.adjust(GO.clust$p, method="BH")

  GO.enriched <- GO.clust[GO.clust$adjustedP <= p, ] %>% arrange(desc(COUNT))

  return(GO.enriched)
}
