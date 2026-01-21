#' Color palette for stacked bar plots
#'
#' A 14-color palette suitable for taxonomic composition plots.
#'
#' @export
stacked_bar.palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
  "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
  "#E4E9B2", "#F9A620", "#054A29", "#52414C"
)

#' Calculate color palette for groups
#'
#' Generates an NPG color palette based on number of groups.
#'
#' @param ps A phyloseq object
#' @param group_variable Name of grouping variable in sample_data
#' @return Character vector of hex colors
#' @export
calc_pal <- function(ps, group_variable) {
  meta <- meta_to_df(ps)
  groups <- unique(meta[, group_variable])

  if (length(groups) < 10) {
    pal <- ggsci::pal_npg()(length(groups))
  } else {
    message("Exceeded colour palette limit, use custom palette")
    pal <- NULL
  }
  return(pal)
}

#' Get top N most abundant taxa
#'
#' @param ps A phyloseq object
#' @param n Number of top taxa to return
#' @param tax_level Taxonomic level (default: "species")
#' @param agg Aggregate at taxonomic level first (default: FALSE)
#' @return Character vector of top taxa names
#' @export
get_top_n <- function(ps, n, tax_level = "species", agg = FALSE) {
  if (tax_level != "species" && agg == TRUE) {
    ps <- ps %>%
      microViz::tax_fix() %>%
      microViz::tax_agg(rank = tax_level)
    ps <- microViz::ps_get(ps)
  }

  topn <- ps %>%
    transform(transform = "relative") %>%
    phyloseq::psmelt() %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise(Mean_abund = mean(Abundance)) %>%
    dplyr::filter(Mean_abund > 0) %>%
    dplyr::slice_max(Mean_abund, n = n) %>%
    dplyr::pull(OTU)

  return(topn)
}

#' Get top N most abundant taxa per group
#'
#' @param ps A phyloseq object
#' @param n Number of top taxa to return
#' @param tax_level Taxonomic level (default: "species")
#' @param var Grouping variable
#' @param group Specific group value to filter
#' @param agg Aggregate at taxonomic level first (default: FALSE)
#' @return Character vector of top taxa names
#' @export
get_top_n_group <- function(ps, n, tax_level = "species", var, group = NULL, agg = FALSE) {
  if (tax_level != "species" && agg != TRUE) {
    warning("Aggregation not set to TRUE, ensure data are aggregated at desired taxonomic level prior to running this function")
  }

  if (tax_level != "species" && agg == TRUE) {
    ps <- ps %>%
      microViz::tax_fix(unknowns = "Unknown") %>%
      microViz::tax_agg(rank = tax_level)
    ps <- microViz::ps_get(ps)
  }

  topn <- ps %>%
    transform(transform = "relative") %>%
    phyloseq::psmelt() %>%
    dplyr::filter({{ var }} == {{ group }}) %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise(Mean_abund = mean(Abundance)) %>%
    dplyr::filter(Mean_abund > 0) %>%
    dplyr::slice_max(Mean_abund, n = n) %>%
    dplyr::pull(OTU)

  return(topn)
}
