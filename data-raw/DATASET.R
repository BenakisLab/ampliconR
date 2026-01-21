## code to prepare `DATASET` dataset goes here
library(ampliconR)
ps_StrokeMice <- import_as_pseq("~/Stroke_Microbiota_reproducibility/data/ASV_seqtab_tax.tab",
                     mapping = "~/Stroke_Microbiota_reproducibility/data/Metadata-16S-sequenced_wo_ctrls.txt",
                     tree="~/Stroke_Microbiota_reproducibility/data/ASV_tree.tre")
usethis::use_data(ps_StrokeMice, overwrite = TRUE)
