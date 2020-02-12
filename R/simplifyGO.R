#' Title
#' @description Remove semantically similar GO terms.
#' @description Based on code by Jamie Soul:
#' @description https://github.com/soulj/SkeletalVis-Pipeline/blob/master/src/SYBIL%20Systems%20Biology/GOEnrichment.R
#'
#' @param GORes dataframe with columns GOID and TERM
#' @param simplifyData output of simplifyGOReqData
#'
#' @return filtered GORes dataframe
#' @importFrom tidyr gather
#' @importFrom GOSemSim mgoSim
#' @importFrom magrittr %>%
#' @export
#'
#' @examples exampleGO <- data.frame(TERM = c(), GOID=c())
#' @examples simplifyGO(exampleGO)
#' 
simplifyGO <- function(GORes, simplifyData = simplifyGOReqData()){
  
  if(nrow(GORes) < 1){
    NA
  }else{
    GORes <- merge(GORes, simplifyData$GOterms, by = "TERM") %>% unique
    
    sim <- mgoSim(GORes$GOID, GORes$GOID,
                  semData = simplifyData$semData,
                  measure="Rel",
                  combine=NULL)
    
    sim[is.na(sim)] <- 0
    
    ## to satisfy codetools for calling gather
    go1 <- go2 <- similarity <- NULL
    
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)
    sim.df <- sim.df[!is.na(sim.df$similarity),]
    sim.df <- sim.df[ order(sim.df$similarity,decreasing = T),]
    
    #get the simiar term pairs
    sim.df <- sim.df[sim.df$go1 !=sim.df$go2,]
    sim.df <- sim.df[sim.df$similarity > 0.4,]
    
    #mark high fequency terms
    sim.df$remove <- apply(sim.df,1,function(x) {
      if (x[1] %in% simplifyData$highFreqTerms){
        return(x[1])
      }
      if (x[2] %in% simplifyData$highFreqTerms){
        return(x[2])
      } else {
        return(NA)
      }
    })
    remove<-na.omit(sim.df$remove)
    sim.df<-sim.df[is.na(sim.df$remove),]
    
    if(nrow(sim.df) == 0){
      NA
    }else{
      #iterate and remove the term with the least significant p-value
      sim.df$go1.pval<-GORes$adjustedP[match(sim.df$go1,GORes$GOID)]
      sim.df$go2.pval<-GORes$adjustedP[match(sim.df$go2,GORes$GOID)]
      
      for (i in 1:nrow(sim.df)){
        
        #check to see if the goterm has already been marked for removal
        if (sim.df[i,"go1"] %in% remove){
          next
        }
        if (sim.df[i,"go2"] %in% remove){
          next
        }
        
        go1.pval <- sim.df[i,"go1.pval"]
        go2.pval <- sim.df[i,"go2.pval"]
        
        #if the p-values are equal then check if parent-child and keep the child term if so
        if (go1.pval==go2.pval){
          go1 <- sim.df[i,"go1"]
          go2 <- sim.df[i,"go2"]
          if(go2 %in% simplifyData$childTerms[[go1]]){
            remove<-c(remove,go2)
            next
          } else if (go1 %in% simplifyData$childTerms[[go2]])
            remove<-c(remove,go1)
          next
        }
        
        #remove least sig term
        remove<-c(remove,sim.df[i,which.max( c(go1.pval,go2.pval))])
      }
      
      GORes.filt <- GORes[ !GORes$GOID %in% remove,]
      
      return(GORes.filt)
    }
  }
}