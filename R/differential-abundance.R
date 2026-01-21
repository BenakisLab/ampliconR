#' Add taxonomy to differential abundance results
#'
#' @param ps A phyloseq object
#' @param list_da_asvs Vector of differentially abundant ASV names
#' @param da_res Differential abundance results data frame
#' @return Data frame with taxonomy added
#' @keywords internal
add_taxonomy_da <- function(ps, list_da_asvs, da_res) {
  taxonomy <- as.data.frame(phyloseq::tax_table(ps))

  da_taxonomy <- taxonomy %>%
    dplyr::select(Highest_classified) %>%
    dplyr::filter(Highest_classified %in% list_da_asvs)

  da_tax_out <- da_res %>%
    tibble::rownames_to_column(var = "Group") %>%
    tibble::column_to_rownames("ASV") %>%
    merge(da_taxonomy, by = 0)

  return(da_tax_out)
}

#' ANCOMBC differential abundance analysis
#'
#' Performs differential abundance testing using ANCOMBC2.
#'
#' @param ps A phyloseq object
#' @param formula Model formula (right-hand side only, e.g., "Group")
#' @param group Grouping variable name
#' @param ord Optional vector specifying group order (first is reference)
#' @param zero_thresh Prevalence threshold for filtering (default: 0.33)
#' @param tax_level Optional taxonomic level for aggregation
#' @param format_tax Whether to format taxonomy before analysis (default: FALSE)
#' @param ... Additional arguments passed to ancombc2()
#' @return Data frame of differentially abundant taxa
#' @export
ancom_da <- function(ps,
                     formula,
                     group,
                     ord = NULL,
                     zero_thresh = 0.33,
                     tax_level = NULL,
                     format_tax = FALSE,
                     ...) {
  if (!is.null(ord)) {
    phyloseq::sample_data(ps)[[group]] <-
      factor(phyloseq::sample_data(ps)[[group]], levels = ord)
  } else {
    ord <- sort(unique(phyloseq::sample_data(ps)[[group]]))
  }

  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

  if (!is.null(tax_level)) {
    if (!tax_level %in% levels) {
      stop(paste0("Not a valid taxonomic rank, please specify from: ",
                  paste(levels, collapse = ", "), ", or leave blank"))
    }
  }

  if (format_tax == TRUE) {
    ps <- format_taxonomy(ps)
  }

  res <- ANCOMBC::ancombc2(
    data = ps,
    tax_level = tax_level,
    fix_formula = formula,
    p_adj_method = "BH",
    prv_cut = zero_thresh,
    group = group,
    struc_zero = TRUE,
    neg_lb = FALSE,
    pseudo_sens = FALSE,
    global = FALSE,
    ...
  )

  res_df <- data.frame(
    Highest_classified = unlist(res$res$taxon),
    lfc = unlist(res$res[3]),
    se = unlist(res$res[5]),
    W = unlist(res$res[7]),
    pval = unlist(res$res[9]),
    qval = unlist(res$res[11]),
    da = unlist(res$res[13])
  ) %>%
    dplyr::mutate(Log2FC = log2(exp(lfc)))

  res_da <- res_df %>%
    dplyr::filter(da == TRUE)

  da_asvs <- res_da$Highest_classified

  if (is.null(tax_level)) {
    if (format_tax == TRUE) {
      res_df <- add_taxonomy_da(ps, da_asvs, res_df)
    }
  }

  reference <- paste(ord[1], "vs", ord[2], sep = "_")
  res_da <- res_da %>%
    dplyr::mutate(Group_ord = reference)

  return(res_da)
}
