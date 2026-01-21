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
