#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gwas_id = args[1]  # gwas gwas_id
pval = as.numeric(args[2])  # gwas p-val
r2 = as.numeric(args[3])  # gwas p-val
kb = as.numeric(args[4])  # gwas p-val
out_tsv_path = args[5]  # output TSV
out_bed_path = args[6]  # output TSV

################################################################################
# Init output df
mat = matrix(ncol = 0, nrow = 0)
out_cols = c("se", "beta", "position", "p", "n", "chr", "id", "rsid", "ea", "nea", "eaf")
out_df = data.frame(matrix(vector(), 0, length(out_cols), dimnames = list(c(), out_cols)), stringsAsFactors = F)


# gwasglue_tophits <- function(gwas_id, pval, out_tsv_path, out_bed_path) {

# Create output dir
dir.create(dirname(out_bed_path), showWarnings = FALSE)

# top gwasglue and convert them to granges
gwas_top_hg19_tbl <- ieugwasr::tophits(gwas_id, pval = pval, clump=1, r2 = r2, kb = kb)
gwas_top_hg19_tbl <- unique(gwas_top_hg19_tbl)
print(head(gwas_top_hg19_tbl))

if (nrow(gwas_top_hg19_tbl)==0) {
    write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
} else if (class(gwas_top_hg19_tbl)[[1]] == "response") { # server code error 503
    write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
} else {
  write.table(gwas_top_hg19_tbl, file = out_tsv_path, sep = "\t", col.names = T, row.names = F, quote = F)
  # Create 4-col bed
  gwas_top_hg19_bed_df = data.frame(chr = gwas_top_hg19_tbl$chr, start = gwas_top_hg19_tbl$position - 1, end = gwas_top_hg19_tbl$position, name = gwas_top_hg19_tbl$rsid)
  write.table(gwas_top_hg19_bed_df, file = out_bed_path, sep = "\t", col.names = F, row.names = F, quote = F)
}
