#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
tophits_hg19_tsv_path = args[1]  #
tophits_hg38_bed_path = args[2]  #
tophits_hg38_tsv_path = args[3]  #

hg19_tsv_df = readr::read_tsv(tophits_hg19_tsv_path, show_col_types = FALSE)
hg38_bed_df = read.table(tophits_hg38_bed_path, col.names=c('chr', 'start', 'position', 'rsid'))

if (nrow(hg19_tsv_df) == 0 | nrow(hg38_bed_df) == 0) { # if either file empty
  cat("Empty input\n")
  cat(NULL, file=tophits_hg38_tsv_path)

} else {

hg38_tsv_df = merge(hg19_tsv_df, hg38_bed_df, by='rsid')
hg38_tsv_df = subset(hg38_tsv_df, select = -c(position.x, chr.x, start))
names(hg38_tsv_df)[names(hg38_tsv_df) == "chr.y"] <- "chr"
names(hg38_tsv_df)[names(hg38_tsv_df) == "position.y"] <- "position"

col_order = c("chr", "position", "rsid", "nea", "ea", "p", "beta", "se", "eaf", "n", "id")
hg38_tsv_df = hg38_tsv_df[, col_order]
hg38_tsv_df = hg38_tsv_df[order(hg38_tsv_df$chr, hg38_tsv_df$position), ]

write.table(hg38_tsv_df, tophits_hg38_tsv_path, quote=F, sep="\t", row.names = F, col.names = T)

}

