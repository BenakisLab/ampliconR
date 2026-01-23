test_that("stacked_bar.palette has 14 colors", {
  expect_length(stacked_bar.palette, 14)
  expect_type(stacked_bar.palette, "character")
})

test_that("stacked_bar.palette contains valid hex colors", {
  # Check that all values are valid hex color codes
  hex_pattern <- "^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$"
  expect_true(all(grepl(hex_pattern, stacked_bar.palette)))
})

test_that("calc_pal() returns color palette for groups", {
  ps <- create_mock_phyloseq()
  pal <- calc_pal(ps, "Group")

  expect_type(pal, "character")
  expect_length(pal, 2)  # Control and Treatment
})

test_that("calc_pal() handles multiple groups", {
  ps <- create_mock_phyloseq()
  pal <- calc_pal(ps, "Timepoint")

  expect_length(pal, 3)  # T1, T2, T3
})

test_that("calc_pal() warns for >9 groups", {
  # Create phyloseq with many groups
  ps <- create_mock_phyloseq()
  # Add 10 unique groups to trigger the warning
  phyloseq::sample_data(ps)$ManyGroups <- paste0("G", c(1:6, 7:10)[1:6])

  # Need to add more samples to have 10 groups
  # Create a larger mock phyloseq with 10 samples
  set.seed(42)
  otu_mat <- matrix(
    sample(0:100, 50, replace = TRUE),
    nrow = 5, ncol = 10,
    dimnames = list(paste0("ASV", 1:5), paste0("Sample", 1:10))
  )
  otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

  sample_df <- data.frame(
    ManyGroups = paste0("G", 1:10),
    row.names = paste0("Sample", 1:10)
  )
  samp <- phyloseq::sample_data(sample_df)

  ps_large <- phyloseq::phyloseq(otu, samp)

  expect_message(pal <- calc_pal(ps_large, "ManyGroups"))
})
