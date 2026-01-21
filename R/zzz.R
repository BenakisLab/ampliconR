#' @import phyloseq
#' @import vegan
#' @import ggplot2
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup pull
#' @importFrom dplyr rowwise c_across across bind_rows distinct slice_max vars
#' @importFrom dplyr row_number if_else case_when mutate_all mutate_at
#' @importFrom tidyr separate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats as.formula anova
#' @importFrom grDevices rgb
#' @importFrom utils read.table
NULL

.onLoad <- function(libname, pkgname) {
}
