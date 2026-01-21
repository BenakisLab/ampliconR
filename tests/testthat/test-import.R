test_that("format_taxonomy() fills missing taxonomy with last known", {
  ps <- create_mock_phyloseq_missing_tax()
  ps_formatted <- format_taxonomy(ps)

  tax_df <- taxonomy(ps_formatted)

  # Check that NAs are filled
  expect_false(any(is.na(tax_df$Genus)))

  # Check that repeated ranks get "unknown_" prefix
  # ASV2 had Clostridia for Class and NA for the rest, should propagate
  asv2_row <- tax_df["ASV2", ]
  expect_true(grepl("unknown_", asv2_row$Order) || asv2_row$Order == asv2_row$Class)
})

test_that("format_taxonomy() creates Highest_classified column", {
  ps <- create_mock_phyloseq()
  ps_formatted <- format_taxonomy(ps)

  tax_df <- taxonomy(ps_formatted)
  expect_true("Highest_classified" %in% colnames(tax_df))
})

test_that("format_taxonomy() updates taxa_names", {
  ps <- create_mock_phyloseq()
  ps_formatted <- format_taxonomy(ps)

  new_names <- phyloseq::taxa_names(ps_formatted)
  # Names should now include genus information
  expect_true(all(grepl(";", new_names) | grepl("ASV", new_names)))
})

test_that("split_taxonomy() splits semicolon-delimited taxonomy", {
  test_df <- data.frame(
    taxonomy = c(
      "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus",
      "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia"
    ),
    row.names = c("ASV1", "ASV2")
  )

  result <- ampliconR:::split_taxonomy(test_df)

  expect_equal(ncol(result), 6)
  expect_equal(colnames(result), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
  expect_equal(result$Kingdom, c("Bacteria", "Bacteria"))
  expect_equal(result$Genus, c("Lactobacillus", "Escherichia"))
})

test_that("read_tab_delim() reads tab-separated files", {
  # Create temporary test file
  tmp_file <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(col1 = 1:3, col2 = c("a", "b", "c"), row.names = c("r1", "r2", "r3")),
    tmp_file,
    sep = "\t",
    quote = FALSE
  )

  result <- ampliconR:::read_tab_delim(tmp_file)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(rownames(result), c("r1", "r2", "r3"))

  unlink(tmp_file)
})

test_that("load_phylo() creates phyloseq from components", {
  # Create components
  otu_mat <- matrix(1:12, nrow = 3, ncol = 4,
                    dimnames = list(paste0("ASV", 1:3), paste0("S", 1:4)))
  tax_mat <- data.frame(
    Kingdom = rep("Bacteria", 3),
    Phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
    row.names = paste0("ASV", 1:3)
  )
  meta <- data.frame(
    Group = c("A", "A", "B", "B"),
    row.names = paste0("S", 1:4)
  )

  result <- ampliconR:::load_phylo(otu_mat, tax_mat, meta)

  expect_s4_class(result, "phyloseq")
  expect_equal(phyloseq::ntaxa(result), 3)
  expect_equal(phyloseq::nsamples(result), 4)
})
