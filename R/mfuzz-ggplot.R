#' Fuzy C-Means Clustering and Visualisation
#' @description  From Jamie Soul; adapted from jspaezp on GitHub
#' @param data a dataframe or matrix
#' @param clustering mFuzz output
#' @param centre logical; should cluster centers be plotted
#' @param sort.columns logical; should column names be replaced by numbered numbers
#'
#' @return list containing ggplot object and dataframe of clustering values
#' @importFrom tidyr gather
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 ggplot aes geom_line facet_grid geom_hline
#' @export
#'
#' @examples data <- matrix(rnorm(1:300), nrow = 100, ncol = 3)
#' @examples 
mFuzz.ggplot <- function(data, clustering,
                         centre = F, sort.columns=T) {

  clusterindex <- clustering$cluster
  
  # data frame with Membership values
  memship <- clustering$membership 
  colnames(memship) <- paste("membership", 
                             seq_along(memship[1,]), 
                             sep = ("")) 
  
  exp <-data
  
  # This chunk replaces col names by numbers if 
  # more than 1 is character only 
  # or when sort.columns is FALSE
  
  all.char.cols <- !grepl("\\d", colnames(exp))
  if ((sum(all.char.cols) > 1) | !sort.columns) {
    colnames(exp) <- seq_along(all.char.cols)    
  }
  
  exp <- data.frame(exp , 
                    Identifier = rownames(data),
                    clusterindex, memship,
                    stringsAsFactors = F) 
  
  # Transform data frame into a ggplot-compatible format
  exp <- exp %>% 
    gather(sample, 
           expression ,
           - Identifier,
           - clusterindex,
           - contains("membership")) %>% 
    mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                       "", 
                       sample,
                       perl = TRUE)) %>%
    #  this regular expression deletes all characters and numbers prior to 
    #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
    mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when seveal character cols ...
                       "0", 
                       Time,
                       perl = TRUE)) %>%
    mutate(Time = factor(Time, ordered=T))
  
  exp$maxMembership <- exp %>%  
    dplyr::select(contains("membership")) %>%
    apply(., 1, max) 
  
  g <- ggplot(exp,
              aes(x = Time, y = expression)) +
    geom_line(
      aes(group = Identifier,  
          colour = maxMembership)
    ) + 
    scale_colour_gradient(low = "white", high = "red")
  
  # Center plotting when centre == TRUE
  if (centre) {
    centers <- clustering$centers %>% 
      data.frame(., 
                 clusterindex = rownames(.),
                 stringsAsFactors = F) %>% 
      gather(sample, 
             Centre,
             - clusterindex) %>% 
      mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                         "", 
                         sample,
                         perl = TRUE)) %>%
      #  this regular expression deletes all characters and numbers prior to 
      #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
      mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when all character names
                         "0", 
                         Time,
                         perl = TRUE)) %>%
      mutate(Time = factor(Time, ordered=T))
    g <- g + geom_line(data = centers, aes(x = Time, y = Centre))
  }
  
  g <- g + facet_grid(. ~ clusterindex)
  g <- g + geom_hline(yintercept = 0, 
                      colour="#666666", alpha = 0.5) + 
    ylab("Fold change") + labs(fill= "Membership")
  
  return(list(g=g,exp=exp))
}
