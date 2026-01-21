#' Extract taxonomy table from phyloseq object
#'
#' @param ps A phyloseq object
#' @return A data frame containing the taxonomy table
#' @export
taxonomy <- function(ps) {

return(as.data.frame(phyloseq::tax_table(ps)))
}

#' Convert sample metadata to data frame
#'
#' @param ps A phyloseq object
#' @return A data frame containing sample metadata
#' @export
meta_to_df <- function(ps) {
  return(as(phyloseq::sample_data(ps), "data.frame"))
}

#' Convert phyloseq OTU table to data frame
#'
#' @param ps A phyloseq object
#' @return A data frame containing the OTU/ASV table
#' @export
ps_to_asvtab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}
