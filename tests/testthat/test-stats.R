test_that("phyloseq_adonis() performs PERMANOVA", {
  ps <- create_mock_phyloseq()
  dist_mat <- phyloseq::distance(ps, "bray")

  result <- phyloseq_adonis(ps, dist_mat, "Group")

  expect_s3_class(result, "anova.cca")
  expect_true("R2" %in% colnames(result))
  expect_true("Pr(>F)" %in% colnames(result))
})

test_that("phyloseq_adonis() returns NULL with NA in grouping variable", {
  ps <- create_mock_phyloseq()
  # Add NA to Group
  phyloseq::sample_data(ps)$Group[1] <- NA
  dist_mat <- phyloseq::distance(ps, "bray")

  expect_message(result <- phyloseq_adonis(ps, dist_mat, "Group"), "NAs")
  expect_null(result)
})

test_that("phyloseq_betadisper() tests homogeneity of dispersion", {
  ps <- create_mock_phyloseq()
  dist_mat <- phyloseq::distance(ps, "bray")

  result <- phyloseq_betadisper(ps, dist_mat, "Group")

  expect_s3_class(result, "anova")
  expect_true("Pr(>F)" %in% colnames(result))
})

test_that("phyloseq_betadisper() returns NULL with NA in grouping variable", {
  ps <- create_mock_phyloseq()
  phyloseq::sample_data(ps)$Group[1] <- NA
  dist_mat <- phyloseq::distance(ps, "bray")

  expect_message(result <- phyloseq_betadisper(ps, dist_mat, "Group"), "NAs")
  expect_null(result)
})

test_that("xyform() creates formula from variable names", {
  form <- xyform("y", "x")
  expect_s3_class(form, "formula")
  expect_equal(deparse(form), "y ~ x")
})

test_that("xyform() handles multiple predictors", {
  form <- xyform("response", c("var1", "var2", "var3"))
  expect_equal(deparse(form), "response ~ var1 + var2 + var3")
})

test_that("round_any() rounds to specified accuracy", {
  expect_equal(round_any(13, 5), 15)
  expect_equal(round_any(13, 5, floor), 10)
  expect_equal(round_any(13, 10), 20)
  expect_equal(round_any(0.123, 0.05), 0.15)
})
