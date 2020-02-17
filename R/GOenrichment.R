#' Over Representation Analysis
#' @description Calculate over enrichment of GO terms in
#' @description a sample Explanation can be found in:
#' @description https://mengnote.blogspot.co.uk/2012/12/calculate-correct-hypergeometric-p.html
#' @param sampleDistr Output of GOdistr() for sample
#' @param sampleSize Size of sample
#' @param bkgdDistr Output of GOdistr() for background population
#' @param bkgdSize Size of background population
#' @param p Significance threshold
#'
#' @return Over represented GO terms as df
#' @importFrom dplyr arrange desc
#' @importFrom stats phyper p.adjust
#' @export
GOenrichment <- function(sampleDistr,
                         sampleSize,
                         bkgdDistr,
                         bkgdSize,
                         p = 0.05){

  successes.sample <- sampleDistr$COUNT    # How many times term 'i' occurred in sample

  successes.bkgd <- bkgdDistr[bkgdDistr$TERM %in% sampleDistr$TERM,] %>%
    arrange(.data$TERM)
  successes.bkgd <- successes.bkgd$COUNT   #How many times term 'i' occured in background

  failure <- bkgdSize - successes.bkgd  #Number of peptides in background minus the peptides with that term

  sampleDistr$p <- phyper(successes.sample-1,
                          successes.bkgd,
                          failure,
                          sampleSize,
                          lower.tail = F)

  sampleDistr$adjustedP <- p.adjust(sampleDistr$p, method="BH")

  GO.enriched <- sampleDistr[sampleDistr$adjustedP <= p, ] %>%
    arrange(desc(.data$COUNT))

  return(GO.enriched)
}
