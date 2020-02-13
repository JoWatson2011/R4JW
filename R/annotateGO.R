#' Map genes to annotated GO terms
#'
#' @param k Vector of gene to map to GOIDs
#'
#' @return Dataframe with GOIDs, terms and genes
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom GO.db GO.db
#' @importFrom AnnotationDbi select
#' @importFrom dplyr select
#' @export
#' @examples annotateGO(c("MAPK1", "MAPK3"))
annotateGO <- function(k){

  GOid <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = k, columns = "GO",
                                keytype = "SYMBOL") %>%
     filter(.data$ONTOLOGY == "BP") %>%
     dplyr::select(.data$SYMBOL, GOID = .data$GO)

  GOterms <- AnnotationDbi::select(GO.db::GO.db,
                                   keys = GOid$GOID,
                                   columns = c("TERM"),
                                   keytype = "GOID")

  annotations <- merge(GOid, GOterms, by = "GOID") %>% unique()

  return(annotations)
}
