#' Calculate Shannon effective diversity
#'
#' @param x Numeric vector of abundances
#' @return Shannon effective (exponential of Shannon entropy)
#' @keywords internal
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <- round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

#' Calculate Faith's phylogenetic diversity
#'
#' @param x Abundance matrix (samples as rows)
#' @param tree Phylogenetic tree
#' @return Faith's PD values
#' @keywords internal
Faiths <- function(x, tree) {
  pd <- picante::pd(x, tree, include.root = FALSE)
  return(pd[, 1])
}

#' Calculate observed richness
#'
#' @param x Numeric vector of abundances
#' @param detection Detection threshold (default: 0.5)
#' @return Number of taxa above detection threshold
#' @export
#' @keywords internal
Richness <- function(x, detection = 0.5) {
  observed <- sum(x > detection)
  return(observed)
}

#' Calculate alpha diversity metrics
#'
#' Computes Richness, Shannon Effective, and Faith's PD for all samples.
#'
#' @param ps A phyloseq object
#' @param ... Additional arguments passed to Richness()
#' @return A data frame with diversity metrics per sample
#' @export
calc_alpha <- function(ps, ...) {
  mat_in <- ps_to_asvtab(ps) %>% t()

  diversity <- setNames(
    data.frame(matrix(ncol = 3, nrow = phyloseq::nsamples(ps))),
    c("Richness", "Shannon.Effective", "Faiths.PD")
  )
  rownames(diversity) <- rownames(meta_to_df(ps))

  diversity$Richness <- apply(mat_in, 1, Richness, ...)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)

  if (is.null(ps@phy_tree)) {
    warning("Phylogenetic tree is required for Faith's PD calculation, returning NA for this metric")
    diversity$Faiths.PD <- NA
  } else {
    tree <- phyloseq::phy_tree(ps)
    diversity$Faiths.PD <- Faiths(mat_in, tree)
  }

  return(diversity)
}
