test_that("taxonomy() extracts taxonomy table as data frame", {
  ps <- create_mock_phyloseq()
  tax_df <- taxonomy(ps)

expect_s3_class(tax_df, "data.frame")
  expect_equal(nrow(tax_df), 5)
  expect_equal(ncol(tax_df), 6)
  expect_true(all(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") %in% colnames(tax_df)))
  expect_equal(rownames(tax_df), paste0("ASV", 1:5))
})

test_that("meta_to_df() extracts sample metadata as data frame", {
  ps <- create_mock_phyloseq()
  meta_df <- meta_to_df(ps)

  expect_s3_class(meta_df, "data.frame")
  expect_equal(nrow(meta_df), 6)
  expect_true(all(c("SampleID", "Group", "Timepoint") %in% colnames(meta_df)))
  expect_equal(rownames(meta_df), paste0("Sample", 1:6))
})

test_that("ps_to_asvtab() extracts OTU table as data frame", {
  ps <- create_mock_phyloseq()
  asv_df <- ps_to_asvtab(ps)

  expect_s3_class(asv_df, "data.frame")
  expect_equal(nrow(asv_df), 5)
  expect_equal(ncol(asv_df), 6)
  expect_equal(rownames(asv_df), paste0("ASV", 1:5))
  expect_equal(colnames(asv_df), paste0("Sample", 1:6))
})

test_that("data wrangling functions preserve data integrity", {
  ps <- create_mock_phyloseq()

  # Check that extracted values match original
  orig_otu <- as.data.frame(phyloseq::otu_table(ps))
  extracted_otu <- ps_to_asvtab(ps)
  expect_equal(orig_otu, extracted_otu)

  orig_tax <- as.data.frame(phyloseq::tax_table(ps))
  extracted_tax <- taxonomy(ps)
  expect_equal(orig_tax, extracted_tax)
})
