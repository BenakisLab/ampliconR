test_that("Shannon.E() calculates Shannon effective diversity", {
  # Test with known values
  # Equal abundances should give highest diversity
  equal_abund <- c(25, 25, 25, 25)
  result <- Shannon.E(equal_abund)

  expect_type(result, "double")
  expect_equal(result, 4.0)  # exp(log(4)) = 4 for equal abundances of 4 species
})

test_that("Shannon.E() handles single species", {
  single <- c(100, 0, 0, 0)
  result <- Shannon.E(single)

  expect_equal(result, 1.0)  # Single species = diversity of 1
})

test_that("Shannon.E() handles zeros correctly", {
  with_zeros <- c(50, 50, 0, 0)
  result <- Shannon.E(with_zeros)

  expect_equal(result, 2.0)  # Two equal species
})

test_that("Richness() counts taxa above threshold", {
  abundances <- c(10, 5, 0.3, 0, 100)

  expect_equal(Richness(abundances, detection = 0.5), 3)
  expect_equal(Richness(abundances, detection = 0), 4)
  expect_equal(Richness(abundances, detection = 50), 1)
})

test_that("Richness() returns 0 for empty samples", {
  empty <- c(0, 0, 0, 0)
  expect_equal(Richness(empty), 0)
})

test_that("Faiths() requires phylogenetic tree", {
  skip_if_not_installed("ape")

  ps <- create_mock_phyloseq(with_tree = TRUE)
  mat <- t(ps_to_asvtab(ps))
  tree <- phyloseq::phy_tree(ps)

  result <- Faiths(mat, tree)

  expect_type(result, "double")
  expect_length(result, 6)  # 6 samples
  expect_true(all(result > 0))
})

test_that("calc_alpha() returns data frame with all metrics", {
  skip_if_not_installed("ape")

  ps <- create_mock_phyloseq(with_tree = TRUE)
  alpha <- calc_alpha(ps)

  expect_s3_class(alpha, "data.frame")
  expect_equal(nrow(alpha), 6)  # 6 samples
  expect_true(all(c("Richness", "Shannon.Effective", "Faiths.PD") %in% colnames(alpha)))
})

test_that("calc_alpha() warns when tree is missing for Faiths PD", {
  ps <- create_mock_phyloseq(with_tree = FALSE)

  expect_warning(alpha <- calc_alpha(ps), "Phylogenetic tree is required")
  expect_true(all(is.na(alpha$Faiths.PD)))
})

test_that("calc_alpha() row names match sample names", {
  ps <- create_mock_phyloseq(with_tree = FALSE)

  suppressWarnings(alpha <- calc_alpha(ps))
  expect_equal(rownames(alpha), paste0("Sample", 1:6))
})
