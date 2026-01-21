#' Calculate Generalized UniFrac distance
#'
#' Wrapper for GUniFrac calculation on phyloseq objects.
#'
#' @param ps A phyloseq object with phylogenetic tree
#' @param asdist Return as dist object (default: TRUE)
#' @return Distance matrix or dist object (alpha = 0.5)
#' @export
phyloseq_gunifrac <- function(ps, asdist = TRUE) {
  asvtab <- t(ps_to_asvtab(ps))
  tree <- phyloseq::phy_tree(ps)

  rooted_tree <- phangorn::midpoint(tree)

  gunifrac <- GUniFrac::GUniFrac(asvtab, rooted_tree,
                                  alpha = c(0.0, 0.5, 1.0))$unifracs

  if (asdist == TRUE) {
    return(stats::as.dist(gunifrac[, , "d_0.5"]))
  } else {
    return(gunifrac[, , "d_0.5"])
  }
}

#' Calculate beta diversity distance and ordination
#'
#' Computes distance matrix and ordination for various beta diversity metrics.
#'
#' @param ps A phyloseq object
#' @param dist Distance metric: "bray", "unifrac", "wunifrac", "gunifrac", or "aitchison"
#' @param ord_method Ordination method: "NMDS", "MDS", or "PCoA" (default: "NMDS")
#' @return A list with Distance_Matrix and Ordination
#' @export
calc_betadiv <- function(ps, dist, ord_method = "NMDS") {
  if (!ord_method %in% c("NMDS", "MDS", "PCoA")) {
    stop("Ordination method not supported. Supported methods: NMDS, MDS, PCoA")
  }

  if (dist %in% c("unifrac", "wunifrac", "gunifrac")) {
    if (is.null(ps@phy_tree)) {
      stop("Phylogenetic tree is required for UniFrac distances")
    }
  }

  if (dist %in% c("unifrac", "wunifrac", "bray")) {
    dist_mat <- phyloseq::distance(ps, dist)
    ord <- phyloseq::ordinate(ps, ord_method, dist_mat)
    return(list("Distance_Matrix" = dist_mat, "Ordination" = ord))
  } else if (dist == "gunifrac") {
    dist_mat <- phyloseq_gunifrac(ps)
    ord <- phyloseq::ordinate(ps, ord_method, dist_mat)
    return(list("Distance_Matrix" = dist_mat, "Ordination" = ord))
  } else if (dist == "aitchison") {
    dist_mat <- phyloseq::distance(transform(ps, "clr"), method = "euclidean")
    ord <- phyloseq::ordinate(ps, ord_method, dist_mat)
    return(list("Distance_Matrix" = dist_mat, "Ordination" = ord))
  } else {
    stop("Distance metric not supported. Supported: bray, unifrac, wunifrac, gunifrac, aitchison")
  }
}
