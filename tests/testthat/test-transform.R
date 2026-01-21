test_that("transform() with mss normalizes to minimum sum", {
  ps <- create_mock_phyloseq()
  ps_mss <- transform(ps, transform = "mss")

  asv_tab <- ps_to_asvtab(ps_mss)
  col_sums <- colSums(asv_tab)

  # All columns should have the same sum (minimum sum scaling)
  expect_true(all(abs(col_sums - min(col_sums)) < 1e-10))
})

test_that("transform() with relative converts to percentages", {
  ps <- create_mock_phyloseq()
  ps_rel <- transform(ps, transform = "relative")

  asv_tab <- ps_to_asvtab(ps_rel)
  col_sums <- colSums(asv_tab)

  # All columns should sum to 100
  expect_true(all(abs(col_sums - 100) < 1e-10))
})

test_that("transform() with clr applies CLR transformation", {
  ps <- create_mock_phyloseq()
  ps_clr <- transform(ps, transform = "clr")

  asv_tab <- ps_to_asvtab(ps_clr)

  # CLR transformed data should have mean of approximately 0 per sample
  col_means <- colMeans(asv_tab)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("transform() returns phyloseq object", {
  ps <- create_mock_phyloseq()

  expect_s4_class(transform(ps, "mss"), "phyloseq")
  expect_s4_class(transform(ps, "relative"), "phyloseq")
  expect_s4_class(transform(ps, "clr"), "phyloseq")
})
test_that("transform() errors on invalid transform type", {
  ps <- create_mock_phyloseq()

  expect_error(transform(ps, transform = "invalid"))
})

test_that("transform() errors on non-phyloseq input", {
  expect_error(transform(data.frame(a = 1:3)))
})
