#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# gwas_id = "ieu-a-834"
# pval = "5e-08"
# r2 = "0.1"
# kb = "1000"
# out_tsv_path = "o.tsv"
# out_bed_path = "o.bed"

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

dir.create(dirname(out_bed_path), showWarnings = FALSE)

gwas_top_hg19_tbl <- ieugwasr::tophits(gwas_id, pval = pval, clump=1, r2 = r2, kb = kb)


write_empty_output <- function() { # create a function with the name my_function
  print(sprintf("gwas_id: %s, nb tophits: %d", gwas_id, 0)) 
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  file.create(out_bed_path)  # empty bed file  
}

if (class(gwas_top_hg19_tbl)[[1]] == "response") {  # 503
  write_empty_output()
} else {
  if (nrow(gwas_top_hg19_tbl)==0) {
    write_empty_output()
  } else {
    gwas_top_hg19_tbl = unique(gwas_top_hg19_tbl)
    print(sprintf("gwas_id: %s, nb tophits: %d", gwas_id, nrow(gwas_top_hg19_tbl)))  
    write.table(gwas_top_hg19_tbl, file = out_tsv_path, sep = "\t", col.names = T, row.names = F, quote = F)
    # Create 4-col bed
    gwas_top_hg19_bed_df = data.frame(chr = gwas_top_hg19_tbl$chr, start = gwas_top_hg19_tbl$position - 1, end = gwas_top_hg19_tbl$position, name = gwas_top_hg19_tbl$rsid)
    write.table(gwas_top_hg19_bed_df, file = out_bed_path, sep = "\t", col.names = F, row.names = F, quote = F)
  }
}
