#' Transform abundance data
#'
#' Normalize phyloseq abundance data using various methods.
#'
#' @param ps A phyloseq object
#' @param transform Transformation method: "mss" (minimum sum scaling),
#'   "relative" (relative abundance), or "clr" (centered log-ratio)
#' @param offset Pseudocount offset for CLR transformation (default: 1)
#' @return A phyloseq object with transformed abundances
#' @export
transform <- function(ps, transform = "mss", offset = 1) {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_asvtab(ps)
  } else {
    stop("Input must be a phyloseq object")
  }

  if (transform %in% c("mss", "relative", "clr")) {
    if (transform == "mss") {
      ps_t <- t(min(colSums(x)) * t(x) / colSums(x))
    } else if (transform == "relative") {
      ps_t <- t(100 * t(x) / colSums(x))
    } else if (transform == "clr") {
      ps_t <- t(mixOmics::logratio.transfo(t(x + offset), logratio = "CLR"))
      class(ps_t) <- "matrix"
    }
    phyloseq::otu_table(ps)@.Data <- ps_t
    return(ps)
  } else {
    stop("Not a valid transform. Supported: mss, relative, clr")
  }
}
