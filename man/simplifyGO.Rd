% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplifyGO.R
\name{simplifyGO}
\alias{simplifyGO}
\title{Title}
\usage{
simplifyGO(GORes, simplifyData = simplifyGOReqData(), simThresh = 0.4)
}
\arguments{
\item{GORes}{dataframe with columns GOID and TERM}

\item{simplifyData}{output of simplifyGOReqData}

\item{simThresh}{threshold for semantic similarity}
}
\value{
filtered GORes dataframe
}
\description{
Remove semantically similar GO terms.

Based on code by Jamie Soul:

https://github.com/soulj/SkeletalVis-Pipeline/blob/master/src/SYBIL%20Systems%20Biology/GOEnrichment.R
}
\examples{
exampleGO <- data.frame(TERM = c(), GOID=c())
simplifyGO(exampleGO)

}
