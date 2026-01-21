#' Plot boxplot with statistical annotations
#'
#' Creates a boxplot with jittered points and statistical significance annotations.
#'
#' @param df Data frame
#' @param variable_col Column name for x-axis grouping variable
#' @param value_col Column name for y-axis values
#' @param comparisons_list List of pairwise comparisons
#' @param fill_var Fill variable (default: variable_col)
#' @param xlab X-axis label (default: variable_col)
#' @param ylab Y-axis label (default: value_col)
#' @param p_title Plot title
#' @param multiple_groups Use Kruskal-Wallis for multiple groups (default: FALSE)
#' @param cols Custom color palette
#' @param group_order Order of groups on x-axis
#' @param paired Use paired tests (default: FALSE)
#' @param ... Additional arguments passed to stat_pvalue_manual()
#' @return A ggplot object
#' @export
plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group_order = NULL,
                         paired = FALSE,
                         ...) {
  if (is.null(cols)) {
    cols <- ggsci::pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")

  if (!is.null(group_order)) {
    df[, variable_col] <- factor(df[, variable_col], levels = group_order)
  }

  formula <- xyform(value_col, variable_col)

  if (multiple_groups == TRUE) {
    if (paired == TRUE) {
      stat_variance <- df %>% rstatix::friedman_test(formula)
      stat_test <- df %>%
        rstatix::pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x = variable_col) %>%
        dplyr::filter(p.adj < 0.05)
    } else {
      stat_variance <- df %>% rstatix::kruskal_test(formula)
      stat_test <- df %>%
        rstatix::pairwise_wilcox_test(formula,
                                      comparisons = comparisons_list,
                                      p.adjust.method = "BH") %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x = variable_col) %>%
        dplyr::filter(p.adj < 0.05)
    }
  } else {
    if (paired == TRUE) {
      stat_test <- df %>%
        rstatix::wilcox_test(formula, paired = TRUE) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x = variable_col) %>%
        dplyr::filter(p < 0.05)
    } else {
      stat_test <- df %>%
        rstatix::wilcox_test(formula) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x = variable_col) %>%
        dplyr::filter(p < 0.05)
    }
  }

  plot <- ggplot2::ggplot(
    df,
    ggplot2::aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    ggplot2::geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    ggplot2::geom_point(size = 1.5,
                        position = ggplot2::position_jitterdodge()) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::stat_boxplot(color = "black",
                          geom = "errorbar",
                          width = 0.2)

  final_plot <- plot +
    ggpubr::theme_classic2() +
    ggplot2::ggtitle(p_title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggpubr::rotate_x_text(angle = 45)

  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  } else {
    datamax <- max(df[, value_col])
    statmax <- max(stat_test["y.position"])

    if (datamax > 100) {
      if (datamax < 250) {
        ybreak <- 50
      } else {
        ybreak <- 100
      }
      ylimit <- round_any(statmax, accuracy = 100)
      ybreakmax <- round_any(datamax, accuracy = 100)
    } else if (datamax < 10) {
      ybreak <- 1
      ylimit <- round_any(statmax, accuracy = 1)
      ybreakmax <- round_any(datamax, accuracy = 1)
    } else {
      ybreak <- 10
      ylimit <- round_any(statmax, accuracy = 10)
      ybreakmax <- round_any(datamax, accuracy = 10)
    }

    if (multiple_groups == TRUE) {
      plot_out <- final_plot +
        ggpubr::stat_pvalue_manual(
          stat_test,
          label = "p.adj.signif",
          hide.ns = TRUE,
          inherit.aes = FALSE,
          ...
        ) +
        ggplot2::scale_y_continuous(breaks = seq(0, ybreakmax, by = ybreak),
                                    limits = c(0, ylimit)) +
        lemon::coord_capped_cart(left = "top", expand = FALSE)
    } else {
      plot_out <- final_plot +
        ggpubr::stat_pvalue_manual(
          stat_test,
          label = "p.signif",
          hide.ns = TRUE,
          inherit.aes = FALSE,
          ...
        ) +
        ggplot2::scale_y_continuous(breaks = seq(0, ybreakmax, by = ybreak),
                                    limits = c(0, ylimit)) +
        lemon::coord_capped_cart(left = "top", expand = FALSE)
    }
  }

  return(plot_out)
}

#' Plot scatter plot with optional correlation
#'
#' @param df Data frame
#' @param x X variable name
#' @param y Y variable name
#' @param point_color Point color
#' @param line_color Regression line color
#' @param fill_color Confidence interval fill color
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param corr_method Correlation method (NULL for no correlation annotation)
#' @param ... Additional arguments passed to stat_cor()
#' @return A ggplot object
#' @export
plot_scatter <- function(df,
                         x,
                         y,
                         point_color,
                         line_color,
                         fill_color,
                         xlab,
                         ylab,
                         corr_method = NULL,
                         ...) {
  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(ggplot2::aes(color = point_color), size = 2.5) +
    ggplot2::geom_smooth(method = "lm",
                         color = line_color,
                         fill = fill_color) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "None",
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 12)
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

  if (!is.null(corr_method)) {
    p <- p + ggpubr::stat_cor(method = corr_method, ...)
    return(p)
  } else {
    return(p)
  }
}

#' Plot beta diversity ordination
#'
#' Creates an ordination plot with PERMANOVA statistics in caption.
#'
#' @param ps A phyloseq object
#' @param dist_matrix Distance matrix
#' @param ordination Ordination object (from calc_betadiv)
#' @param group_variable Grouping variable name
#' @param add_ellipse Add confidence ellipses (default: FALSE)
#' @param cols Custom color palette
#' @return A ggplot object
#' @export
plot_beta_div <- function(ps,
                          dist_matrix,
                          ordination,
                          group_variable,
                          add_ellipse = FALSE,
                          cols = NULL) {
  if (is.null(cols)) {
    cols <- calc_pal(ps, group_variable)
  }

  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  betadisp <- phyloseq_betadisper(ps, dist_matrix, group_variable)

  if (dim(phyloseq::sample_data(ps))[2] < 2) {
    phyloseq::sample_data(ps)[, 2] <- phyloseq::sample_data(ps)[, 1]
  }

  plot <- phyloseq::plot_ordination(ps, ordination, color = group_variable)
  plot$layers[[1]] <- NULL

  plot_out <- plot +
    ggplot2::geom_point(size = 3, alpha = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(caption = bquote(Adonis ~ R^2 ~ .(round(ad$R2[1], 2)) ~
                                     ~ p - value ~ .(ad$`Pr(>F)`[1])))

  if (add_ellipse == TRUE) {
    plot_out <- plot_out +
      ggplot2::geom_polygon(stat = "ellipse",
                            ggplot2::aes(fill = .data[[group_variable]]),
                            alpha = 0.3)
  }

  if (betadisp$`Pr(>F)`[[1]] < 0.05) {
    warning("Group dispersion is not homogenous, interpret results carefully",
            call. = FALSE)
  }

  return(plot_out)
}

#' Plot differential abundance results
#'
#' Creates a lollipop plot showing log2 fold changes of DA taxa.
#'
#' @param ancom_res Results from ancom_da()
#' @param groups Character vector of group names (length 2)
#' @param cols Custom color palette (length 2)
#' @return A ggplot object (also prints the plot)
#' @export
plot_da <- function(ancom_res, groups, cols = NULL) {
  if (is.null(cols)) {
    cols <- ggsci::pal_npg()(length(groups))
  }

  ancom_da_plot <- ancom_res %>%
    dplyr::mutate(enriched_in = ifelse(Log2FC > 0, groups[2], groups[1]))

  ancom_da_plot_sort <- ancom_da_plot %>%
    dplyr::arrange(Log2FC) %>%
    dplyr::mutate(Highest_classified = factor(Highest_classified, levels = Highest_classified))

  p <- ggplot2::ggplot(
    ancom_da_plot_sort,
    ggplot2::aes(x = Highest_classified, y = Log2FC, color = enriched_in)
  ) +
    ggplot2::geom_point(size = 5) +
    ggplot2::labs(y = paste("\nLog2 Fold-Change", groups[1], "vs", groups[2], sep = " "),
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", size = 14),
      axis.text.y = ggplot2::element_text(color = "black", size = 14),
      axis.title.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14),
      legend.position = "none"
    ) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::scale_color_manual(values = cols)

  print(p)
}


#' Plot taxonomic composition
#'
#' Creates stacked bar plot of taxonomic composition.
#'
#' @param ps A phyloseq object
#' @param tax_level Taxonomic level for plotting
#' @param var Faceting variable
#' @param ord Order of facet levels
#' @param n_taxa Number of top taxa to show (default: 10)
#' @param per_group Get top N per group instead of overall (default: FALSE)
#' @param groups Groups for per-group top N selection
#' @param agg Aggregate at taxonomic level (default: FALSE)
#' @return A ggplot object
#' @export
plot_taxonomic_comp <- function(ps,
                                tax_level,
                                var,
                                ord = NULL,
                                n_taxa = 10,
                                per_group = FALSE,
                                groups = NULL,
                                agg = FALSE) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  if (!tax_level %in% ranks) {
    stop("Provide a valid taxonomic rank")
  }

  var <- rlang::enquo(var)

  if (per_group == TRUE) {
    if (is.null(groups)) {
      groups <- unique(meta_to_df(ps) %>% dplyr::pull(.data[[var]]))
    }
    top_tax <- purrr::map(groups, function(x) {
      ps %>%
        get_top_n_group(
          n = n_taxa,
          level = tax_level,
          var = .data[[var]],
          group = x,
          agg = agg
        )
    })
    top_tax <- c(unique(unlist(top_tax)), "other")
    n_taxa <- length(top_tax)
  } else {
    top_tax <- c(get_top_n(ps, n = n_taxa, level = tax_level), "other")
  }

  if (!is.null(ord)) {
    ps <- ps %>%
      phyloseq::ps_mutate(plot_var = factor(.data[[var]], levels = ord))
  }

  if (agg == TRUE) {
    melt_ps <- ps %>%
      microViz::aggregate_taxa(level = tax_level) %>%
      phyloseq::psmelt()
  } else {
    melt_ps <- ps %>%
      phyloseq::psmelt()
  }

  taxa_plot <- melt_ps %>% dplyr::filter(OTU %in% top_tax)

  ranks <- ranks[ranks %in% colnames(melt_ps)]
  repl_cols <- c("OTU", ranks)

  other_plot <- melt_ps %>%
    dplyr::filter(!OTU %in% top_tax) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Abundance = sum(Abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Sample, .keep_all = TRUE) %>%
    dplyr::mutate_at(dplyr::vars(repl_cols), ~ paste("other"))

  plot_df <- dplyr::bind_rows(taxa_plot, other_plot)

  taxa_pal <- c(stacked_bar.palette[1:(length(top_tax) - 1)], "#DCDCDC")
  names(taxa_pal) <- plot_df %>%
    dplyr::filter(OTU %in% top_tax) %>%
    dplyr::pull(.data[[tax_level]]) %>%
    unique()

  comp_fig <- plot_df %>%
    ggplot2::ggplot(ggplot2::aes(fill = .data[[tax_level]], y = Abundance, x = Sample)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_fill_manual(values = taxa_pal) +
    ggplot2::facet_grid(stats::as.formula(paste("~", var)),
                        scales = "free",
                        space = "free") +
    cowplot::theme_cowplot() +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 22),
      axis.text.y = ggplot2::element_text(size = 18),
      axis.ticks = ggplot2::element_blank()
    )

  return(comp_fig)
}
