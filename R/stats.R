#' PERMANOVA on phyloseq object
#'
#' Runs adonis2 (PERMANOVA) on a phyloseq object with specified distance matrix.
#'
#' @param ps A phyloseq object
#' @param dist_matrix Distance matrix (from calc_betadiv or distance())
#' @param group_variable Name of grouping variable in sample_data
#' @param ... Additional arguments passed to adonis2()
#' @return adonis2 result object
#' @export
phyloseq_adonis <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    message("metadata contains NAs, remove these samples with subset_samples before continuing")
    return(NULL)
  } else {
    form <- as.formula(paste("dist_matrix", group_variable, sep = " ~ "))
    ps_ad <- vegan::adonis2(form, data = meta_df, ...)
    return(ps_ad)
  }
}

#' Test for homogeneity of group dispersions
#'
#' Runs betadisper and ANOVA to test for differences in group dispersion.
#'
#' @param ps A phyloseq object
#' @param dist_matrix Distance matrix
#' @param group_variable Name of grouping variable in sample_data
#' @param ... Additional arguments passed to betadisper()
#' @return ANOVA result for betadisper
#' @export
phyloseq_betadisper <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    message("metadata contains NAs, remove these samples with subset_samples before continuing")
    return(NULL)
  } else {
    bd <- vegan::betadisper(dist_matrix, meta_df[[group_variable]])
    anova_res <- stats::anova(bd)
    return(anova_res)
  }
}

#' Create formula from variable names
#'
#' @param y_var Response variable name
#' @param x_vars Character vector of predictor variable names
#' @return A formula object
#' @keywords internal
xyform <- function(y_var, x_vars) {
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

#' Round to nearest value
#'
#' @param x Numeric value to round
#' @param accuracy Rounding accuracy
#' @param f Rounding function (default: ceiling)
#' @return Rounded value
#' @keywords internal
round_any <- function(x, accuracy, f = ceiling) {
  f(x / accuracy) * accuracy
}
