# Helper file to create mock phyloseq objects for testing

library(phyloseq)

#' Create a minimal mock phyloseq object for testing
#' @param with_tree Include a phylogenetic tree (default: FALSE)
#' @return A phyloseq object
create_mock_phyloseq <- function(with_tree = FALSE) {
  # Create OTU table (5 ASVs x 6 samples)
  set.seed(42)
  otu_mat <- matrix(
    sample(0:100, 30, replace = TRUE),
    nrow = 5, ncol = 6,
    dimnames = list(
      paste0("ASV", 1:5),
      paste0("Sample", 1:6)
    )
  )
  otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

  # Create taxonomy table
  tax_mat <- matrix(
    c(
      "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus",
      "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Blautia",
      "Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides",
      "Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia",
      "Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium"
    ),
    nrow = 5, ncol = 6, byrow = TRUE,
    dimnames = list(
      paste0("ASV", 1:5),
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    )
  )
  tax <- tax_table(tax_mat)

  # Create sample metadata
  sample_df <- data.frame(
    SampleID = paste0("Sample", 1:6),
    Group = rep(c("Control", "Treatment"), each = 3),
    Timepoint = rep(c("T1", "T2", "T3"), 2),
    row.names = paste0("Sample", 1:6)
  )
  samp <- sample_data(sample_df)

  if (with_tree) {
    # Create a simple random tree with meaningful branch lengths
    tree <- ape::rtree(5, tip.label = paste0("ASV", 1:5), br = stats::runif)
    ps <- phyloseq(otu, tax, samp, tree)
  } else {
    ps <- phyloseq(otu, tax, samp)
  }

  return(ps)
}

#' Create mock phyloseq with missing taxonomy
#' @return A phyloseq object with some NA taxonomy values
create_mock_phyloseq_missing_tax <- function() {
  set.seed(42)
  otu_mat <- matrix(
    sample(0:100, 15, replace = TRUE),
    nrow = 3, ncol = 5,
    dimnames = list(
      paste0("ASV", 1:3),
      paste0("Sample", 1:5)
    )
  )
  otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

  # Taxonomy with missing values
  tax_mat <- matrix(
    c(
      "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus",
      "Bacteria", "Firmicutes", "Clostridia", NA, NA, NA,
      "Bacteria", NA, NA, NA, NA, NA
    ),
    nrow = 3, ncol = 6, byrow = TRUE,
    dimnames = list(
      paste0("ASV", 1:3),
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    )
  )
  tax <- tax_table(tax_mat)

  sample_df <- data.frame(
    SampleID = paste0("Sample", 1:5),
    Group = c("A", "A", "B", "B", "B"),
    row.names = paste0("Sample", 1:5)
  )
  samp <- sample_data(sample_df)

  phyloseq(otu, tax, samp)
}

#' Create mock data frame for plotting tests
#' @return A data frame suitable for plot_boxplot tests
create_mock_plot_df <- function() {
  set.seed(42)
  data.frame(
    Group = rep(c("Control", "Treatment"), each = 10),
    Value = c(rnorm(10, mean = 5, sd = 1), rnorm(10, mean = 7, sd = 1)),
    stringsAsFactors = FALSE
  )
}
