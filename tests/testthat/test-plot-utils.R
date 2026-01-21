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
  # Add more unique timepoints
  phyloseq::sample_data(ps)$ManyGroups <- paste0("G", 1:6)

  # Duplicate samples to get more groups
  expect_message(pal <- calc_pal(ps, "ManyGroups"))
})
