#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
identifier = args[1]  # gwas identifier
pval = as.numeric(args[2])  # gwas p-val
r2 = as.numeric(args[3])  # gwas p-val
kb = as.numeric(args[4])  # gwas p-val
out_tsv_path = args[5]  # output TSV
out_bed_path = args[6]  # output BEd

# gwasglue_tophits <- function(identifier, pval, out_tsv_path, out_bed_path) {

# Create output dir
dir.create(dirname(out_bed_path), showWarnings = FALSE)

# top gwasglue and convert them to granges
gwas_top_hg19_tbl <- ieugwasr::tophits(identifier, pval = pval, r2 = r2, kb = kb)
if (class(gwas_top_hg19_tbl)[[1]] == "response") { # server code error 503
  cat("No significant GWAS signals at p-value: ", sprintf("%f", pval), "\n")
  cat(NULL, file = out_tsv_path)
  cat(NULL, file = out_bed_path)
} else if (class(gwas_top_hg19_tbl)[[1]] == "tbl_df" && nrow(gwas_top_hg19_tbl) == 0) { # if non top gwas, then write empty file
  cat("No significant GWAS signals at p-value: ", sprintf("%f", pval), "\n")
  cat(NULL, file = out_tsv_path)
  cat(NULL, file = out_bed_path)
} else {
  write.table(gwas_top_hg19_tbl, file = out_tsv_path, sep = "\t", col.names = T, row.names = F, quote = F)
  # Create 4-col bed
  gwas_top_hg19_bed_df = data.frame(chr = gwas_top_hg19_tbl$chr, start = gwas_top_hg19_tbl$position - 1, end = gwas_top_hg19_tbl$position, name = gwas_top_hg19_tbl$rsid)
  write.table(gwas_top_hg19_bed_df, file = out_bed_path, sep = "\t", col.names = F, row.names = F, quote = F)
}
# }

# gwasglue_tophits(identifier, pval, out_tsv_path, out_bed_path)

