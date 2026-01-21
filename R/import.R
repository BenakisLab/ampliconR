#' Split semicolon-delimited taxonomy string into rank columns
#'
#' @param asvtab A data frame with a 'taxonomy' column
#' @return A data frame with taxonomy split into Kingdom through Genus columns
#' @keywords internal
split_taxonomy <- function(asvtab) {
  taxonomy_phylo <- dplyr::select(asvtab, "taxonomy") %>%
    tidyr::separate("taxonomy",
                    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                    sep = ";")
  return(taxonomy_phylo)
}

#' Format and clean taxonomy table
#'
#' Fills missing taxonomic ranks with last known classification and adds
#' 'unknown_' prefix where ranks match parent rank. Creates Highest_classified
#' column combining ASV ID with lowest classified rank.
#'
#' @param ps A phyloseq object
#' @return A phyloseq object with formatted taxonomy
#' @export
format_taxonomy <- function(ps) {
  ps_tax <- taxonomy(ps) %>%
    dplyr::mutate_all(dplyr::na_if, "")

  ps_tax[] <- t(zoo::na.locf(t(ps_tax))) %>%
    as.data.frame()

  ps_tax <- ps_tax %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dplyr::across(Kingdom:Genus,
                                ~ dplyr::if_else(all(is.na(dplyr::c_across(Kingdom:Genus))),
                                                 "Unclassified",
                                                 .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Genus = dplyr::case_when(Genus == Family ~ paste0("unknown_", Genus),
                                           TRUE ~ Genus)) %>%
    dplyr::mutate(Family = dplyr::case_when(Family == Order ~ paste0("unknown_", Family),
                                            TRUE ~ Family)) %>%
    dplyr::mutate(Order = dplyr::case_when(Order == Class ~ paste0("unknown_", Order),
                                           TRUE ~ Order)) %>%
    dplyr::mutate(Class = dplyr::case_when(Class == Phylum ~ paste0("unknown_", Class),
                                           TRUE ~ Class)) %>%
    dplyr::mutate(Phylum = dplyr::case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                                            TRUE ~ Phylum)) %>%
    dplyr::mutate(ASV = rownames(ps_tax)) %>%
    dplyr::mutate(Highest_classified = dplyr::case_when(
      Genus == "Unclassified" ~ rownames(ps_tax),
      TRUE ~ paste(rownames(ps_tax), Genus, sep = "; ")
    )) %>%
    as.data.frame()

  rownames(ps_tax) <- ps_tax[, "Highest_classified"]
  phyloseq::taxa_names(ps) <- ps_tax[, "Highest_classified"]
  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(ps_tax))

  return(ps)
}

#' Read tab-delimited file
#'
#' @param df Path to tab-delimited file
#' @return A data frame
#' @keywords internal
read_tab_delim <- function(df) {
  df_out <- read.table(
    df,
    row.names = 1,
    header = 1,
    sep = "\t",
    comment.char = "",
    check.names = FALSE
  )
  return(df_out)
}

#' Load data into phyloseq object
#'
#' @param asvtab ASV/OTU abundance matrix
#' @param taxa Taxonomy table
#' @param mapping Sample metadata
#' @param tree Optional phylogenetic tree file path
#' @return A phyloseq object
#' @keywords internal
load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  phylo_asv <- phyloseq::otu_table(asvtab, taxa_are_rows = TRUE)
  phylo_tax <- phyloseq::tax_table(as.matrix(taxa))
  phylo_map <- phyloseq::sample_data(mapping)

  if (!is.null(tree)) {
    phylo_tree <- phyloseq::read_tree(tree)
    return(phyloseq::merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  } else {
    return(phyloseq::merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}

#' Import ASV table, metadata, and optional tree into phyloseq
#'
#' Main import function for creating phyloseq objects from standard
#' tab-delimited files.
#'
#' @param asvtab Path to ASV table (tab-delimited, with 'taxonomy' column)
#' @param mapping Path to sample metadata (tab-delimited)
#' @param tree Optional path to phylogenetic tree file
#' @return A phyloseq object
#' @export
import_as_pseq <- function(asvtab, mapping, tree = NULL) {
  asvtab_taxa <- read_tab_delim(asvtab)
  metadata <- read_tab_delim(mapping)

  asv_matrix <- as.matrix(subset(asvtab_taxa, select = -taxonomy))
  taxonomy_phylo <- split_taxonomy(asvtab_taxa)
  rownames(taxonomy_phylo) <- rownames(asvtab_taxa)

  n_unclassified <- taxonomy_phylo %>%
    dplyr::filter(dplyr::if_all(Kingdom:Genus, is.na)) %>%
    nrow()

  if (n_unclassified > 0) {
    warning(sprintf("%d ASV(s) have no taxonomic classification at any rank, run format_taxonomy to fix and ensure this is expected", n_unclassified))
  }

  if (!is.null(tree)) {
    out <- load_phylo(asv_matrix, taxonomy_phylo, metadata, tree = tree)
  } else {
    out <- load_phylo(asv_matrix, taxonomy_phylo, metadata)
  }
  return(out)
}
