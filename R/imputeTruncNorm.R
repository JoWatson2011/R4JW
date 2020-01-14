#' Title Impute Data by drawing randomly from a distribution
#' @description The default represents a stimulation to represent a population
#' @describeton below the limits of detection, as described in Tyanova et al., 2016
#' @param v a vector containing NA
#'
#' @return imputed values
#' @importFrom msm rtnorm
#' @importFrom stats sd
#' @export
#'
#' @examples foo <- rnorm(100)
#' @examples foo[sample(100, 5)] <- NA
#' @examples imputeTruncNorm(foo)
#'
imputeTruncNorm <- function(v){
  dis <- rtnorm(length(v),
                lower = -(sd(v, na.rm = T)),
                upper = (-(sd(v, na.rm = T)) + 1.0)
  )
  v[is.na(v)] <- sample(dis, length(v[is.na(v)]), replace = F)
  return(v)
}
