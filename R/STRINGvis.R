#' Construct a ggraph visualisation and igraph object
#' from HGNC names and STRING network.
#'
#' @param genes a vector of HGNC identifiers
#' @param stringNW Output of importSTRINGnw
#'
#' @return list of ggraph visualisation and igraph network
#' @importFrom igraph graph_from_data_frame simplify vcount
#' @importFrom ggraph ggraph geom_edge_link theme_graph geom_node_label geom_node_point
#' @importFrom ggplot2 aes
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @export
#'
#' @examples nw <- importSTRINGnw(400)
#' @examples STRINGvis(c("MAPK1","MAPK3"), nw)
STRINGvis <- function(genes, stringNW){
  pal3 <- brewer.pal(10, "Set3")
  set.seed(1)
  nodes <- stringNW[(stringNW$protein1 %in% genes) &
                      (stringNW$protein2 %in% genes), c(2:3)] %>%
    na.omit
  nw <- graph_from_data_frame(nodes, directed = F) %>%
    simplify(remove.multiple = F,remove.loops = T)
  if (vcount(nw) == 0) {
    return (NA)
  }
  g <- suppressMessages(ggraph(nw, layout = "kk") +
                          geom_edge_link(colour = pal3[1], edge_width = 1) +
                          geom_node_point(size = 4) +
                          theme_graph() +
                          geom_node_label(aes(label = name), size = 2.5)
  )
  return(list(g=g, nw=nw))
}
