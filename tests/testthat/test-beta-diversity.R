test_that("calc_betadiv() returns list with distance matrix and ordination", {
  ps <- create_mock_phyloseq()
  result <- calc_betadiv(ps, dist = "bray", ord_method = "NMDS")

  expect_type(result, "list")
  expect_named(result, c("Distance_Matrix", "Ordination"))
  expect_s3_class(result$Distance_Matrix, "dist")
})

test_that("calc_betadiv() works with different distance metrics", {
  ps <- create_mock_phyloseq()

  # Bray-Curtis
  result_bray <- calc_betadiv(ps, dist = "bray")
  expect_s3_class(result_bray$Distance_Matrix, "dist")

  # Aitchison
  result_ait <- calc_betadiv(ps, dist = "aitchison")
  expect_s3_class(result_ait$Distance_Matrix, "dist")
})

test_that("calc_betadiv() works with different ordination methods", {
  ps <- create_mock_phyloseq()

  result_nmds <- calc_betadiv(ps, dist = "bray", ord_method = "NMDS")
  result_pcoa <- calc_betadiv(ps, dist = "bray", ord_method = "PCoA")

  expect_type(result_nmds$Ordination, "list")
  expect_type(result_pcoa$Ordination, "list")
})

test_that("calc_betadiv() errors on invalid distance metric", {
  ps <- create_mock_phyloseq()

  expect_error(calc_betadiv(ps, dist = "invalid"))
})

test_that("calc_betadiv() errors on invalid ordination method", {
  ps <- create_mock_phyloseq()

  expect_error(calc_betadiv(ps, dist = "bray", ord_method = "invalid"))
})

test_that("calc_betadiv() requires tree for UniFrac distances", {
  ps <- create_mock_phyloseq(with_tree = FALSE)

  expect_error(calc_betadiv(ps, dist = "unifrac"), "Phylogenetic tree is required")
  expect_error(calc_betadiv(ps, dist = "wunifrac"), "Phylogenetic tree is required")
  expect_error(calc_betadiv(ps, dist = "gunifrac"), "Phylogenetic tree is required")
})

test_that("calc_betadiv() works with UniFrac when tree present", {
  skip_if_not_installed("ape")

  ps <- create_mock_phyloseq(with_tree = TRUE)

  result_unifrac <- calc_betadiv(ps, dist = "unifrac")
  expect_s3_class(result_unifrac$Distance_Matrix, "dist")

  result_wunifrac <- calc_betadiv(ps, dist = "wunifrac")
  expect_s3_class(result_wunifrac$Distance_Matrix, "dist")
})

test_that("phyloseq_gunifrac() returns distance object", {
  skip_if_not_installed("ape")

  ps <- create_mock_phyloseq(with_tree = TRUE)
  result <- phyloseq_gunifrac(ps, asdist = TRUE)

  expect_s3_class(result, "dist")
  expect_equal(attr(result, "Size"), 6)  # 6 samples
})

test_that("phyloseq_gunifrac() can return matrix", {
  skip_if_not_installed("ape")

  ps <- create_mock_phyloseq(with_tree = TRUE)
  result <- phyloseq_gunifrac(ps, asdist = FALSE)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(6, 6))
})
