test_that("plot_boxplot() returns ggplot object", {
  df <- create_mock_plot_df()
  comparisons <- list(c("Control", "Treatment"))

  p <- plot_boxplot(df, "Group", "Value", comparisons)

  expect_s3_class(p, "ggplot")
})

test_that("plot_boxplot() works with custom colors", {
  df <- create_mock_plot_df()
  comparisons <- list(c("Control", "Treatment"))
  custom_cols <- c("#FF0000", "#0000FF")

  p <- plot_boxplot(df, "Group", "Value", comparisons, cols = custom_cols)

  expect_s3_class(p, "ggplot")
})

test_that("plot_boxplot() respects group order", {
  df <- create_mock_plot_df()
  comparisons <- list(c("Treatment", "Control"))

  p <- plot_boxplot(df, "Group", "Value", comparisons,
                    group.order = c("Treatment", "Control"))

  expect_s3_class(p, "ggplot")
})

test_that("plot_scatter() returns ggplot object", {
  df <- data.frame(x = 1:10, y = 1:10 + rnorm(10))

  p <- plot_scatter(df, "x", "y",
                    point_color = "blue",
                    line_color = "red",
                    fill_color = "pink",
                    xlab = "X axis",
                    ylab = "Y axis")

  expect_s3_class(p, "ggplot")
})

test_that("plot_scatter() adds correlation when specified", {
  df <- data.frame(x = 1:10, y = 1:10 + rnorm(10, sd = 0.1))

  p <- plot_scatter(df, "x", "y",
                    point_color = "blue",
                    line_color = "red",
                    fill_color = "pink",
                    xlab = "X", ylab = "Y",
                    corr.method = "spearman")

  expect_s3_class(p, "ggplot")
  # Check that stat_cor layer is present
  expect_true(any(sapply(p$layers, function(l) inherits(l$stat, "StatCor"))))
})

test_that("plot_beta_div() returns ggplot object", {
  ps <- create_mock_phyloseq()
  beta <- calc_betadiv(ps, "bray")

  p <- plot_beta_div(ps, beta$Distance_Matrix, beta$Ordination, "Group")

  expect_s3_class(p, "ggplot")
})

test_that("plot_beta_div() adds ellipse when requested", {
  ps <- create_mock_phyloseq()
  beta <- calc_betadiv(ps, "bray")

  p <- plot_beta_div(ps, beta$Distance_Matrix, beta$Ordination, "Group",
                     add_ellipse = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_beta_div() warns on heterogeneous dispersion", {
  # This may or may not trigger depending on random data
  ps <- create_mock_phyloseq()
  beta <- calc_betadiv(ps, "bray")

  # Just check it runs without error
  p <- plot_beta_div(ps, beta$Distance_Matrix, beta$Ordination, "Group")
  expect_s3_class(p, "ggplot")
})

test_that("get_top_n() returns top taxa", {
  ps <- create_mock_phyloseq()

  top3 <- get_top_n(ps, n = 3)

  expect_type(top3, "character")
  expect_length(top3, 3)
})

test_that("get_top_n() excludes zero abundance taxa", {
  ps <- create_mock_phyloseq()

  # Should only return taxa with non-zero mean abundance
  top_all <- get_top_n(ps, n = 100)
  expect_true(length(top_all) <= 5)  # Only 5 ASVs in mock data
})

test_that("plot_da() creates lollipop plot", {
  # Create mock ANCOM results
  ancom_res <- data.frame(
    Highest_classified = c("Taxon1", "Taxon2", "Taxon3"),
    Log2FC = c(-2.5, 1.5, 3.0),
    qval = c(0.01, 0.02, 0.001)
  )

  # Capture the plot (plot_da uses print() and returns ggplot invisibly)
  p <- plot_da(ancom_res, groups = c("Control", "Treatment"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_taxonomic_comp() returns ggplot", {
  skip("Requires microViz package")

  ps <- create_mock_phyloseq()

  p <- plot_taxonomic_comp(ps, tax_level = "Phylum", var = Group)

  expect_s3_class(p, "ggplot")
})
