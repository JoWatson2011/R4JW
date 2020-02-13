#' Title Create network from STRING databases
#'
#' @param conf Confidence threshold for STRING edges
#' @param speciesID NCBI species ID. Default Homo sapiens, 9606
#'
#' @return Edgelist of human STRING network
#' @importFrom biomaRt useMart getBM
#' @importFrom STRINGdb STRINGdb
#' @importFrom dplyr mutate select filter mutate_all
#' @importFrom igraph as_data_frame
#' @importFrom rlang .data
#' @export
importSTRINGnw <- function(conf, speciesID = 9606){
  string_db <- STRINGdb$new(version = "10", species = speciesID, score_threshold = 400, input_directory = "")
  tmp <- string_db$load_all()
  STRING.expmt <- as_data_frame(tmp, what = "edges") %>%
    filter(.data$experiments > conf) %>%
    select(protein1 = .data$from, protein2 = .data$to) %>%
    mutate_all(~ sub("9606.", "", .data, fixed = T)) %>%
    mutate(id = 1:nrow(.data))

  #ENSEMBL --> GENE name
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl" )
  genes1 <- unique(c(STRING.expmt$protein1, STRING.expmt$protein2))
  G_list <- getBM(filters = "ensembl_peptide_id",
                  attributes= c("ensembl_peptide_id", "hgnc_symbol"),
                  values=genes1,
                  mart= ensembl)

  #match ensembl ids to gene symbols
  STRING.prot1 <- merge(STRING.expmt[,c("protein1", "id")],
                        G_list,
                        by.x="protein1",
                        by.y ="ensembl_peptide_id",
                        sort = T)
  STRING.prot2 <- merge(STRING.expmt[,c("protein2", "id")],
                        G_list,
                        by.x="protein2",
                        by.y ="ensembl_peptide_id",
                        sort = F)
  STRING.expmt.gene <- merge(STRING.prot1,
                             STRING.prot2,
                             by = "id",
                             all = F) %>%
    filter(!duplicated(.data$id))

  #remove ensembl ids
  STRING.expmt.gene <- STRING.expmt.gene[,c(2,4)]
  colnames(STRING.expmt.gene) <- c("protein1","protein2")

  return(STRING.expmt.gene)
}
