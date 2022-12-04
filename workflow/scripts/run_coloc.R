suppressWarnings(suppressPackageStartupMessages({
  library(coloc)
  library(DBI)
  library(dplyr)
  library(gwasvcf)
  library(VariantAnnotation)
}))

# PARAMS
gwas_id = "ieu-a-801"
eqtl_id = "Kasela_2017_microarray_T-cell_CD8"
window = 1000000
pval = 5e-08
tophits_tsv_path = "out/gwas418/tophits/ieu-a-801/pval_5e-08/r2_0.1/kb_1000/hg38.tsv"
gwas_vcf_path = "/home/gonzalez/Software/process/hg38/gwas.mrcieu.ac.uk/files/ieu-a-801/ieu-a-801.vcf.bgz"
eqtl_permuted_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.permuted.tsv.gz"
eqtl_all_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.all.tsv.gz"
eur_af_sqlite = "out/eur_af.sqlite"
out_tsv_path = "out/gwas418/coloc/ieu-a-801/pval_5e-08/r2_0.1/kb_1000/window_1000000/Kasela_2017_microarray_T-cell_CD8.tsv"
# END PARAMS

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=10) {
  stop("10 arguments must be supplied", call.=FALSE)
}
gwas_id = args[1]
eqtl_id = args[2]
window = as.numeric(args[3])
pval = as.numeric(args[4])
tophits_tsv_path = args[5]
gwas_vcf_path = args[6]
eqtl_permuted_path = args[7]
eqtl_all_path = args[8]
eur_af_sqlite = args[9]
out_tsv_path = args[10]

dir.create(dirname(out_tsv_path), showWarnings = FALSE, recursive = TRUE)

# no bcftools error message and exit
if (length(which(!is.na(system('which bcftools', intern = T)))) == 0) {
  stop("Error. Please, install bcftools")
} else {
  bcftools_path = system("which bcftools", intern = T)
  gwasvcf::set_bcftools(bcftools_path)
}

################################################################################
# Init output df
mat = matrix(ncol = 0, nrow = 0)
out_cols = c(
  "chrom",
  "pos",
  "rsid",
  "ref",
  "alt",
  "eqtl_gene_id",
  "gwas_beta",
  "gwas_pval",
  "gwas_id",
  "eqtl_beta",
  "eqtl_pval",
  "eqtl_id",
  "PP.H4.abf",
  "SNP.PP.H4",
  "nsnps",
  "PP.H3.abf",
  "PP.H2.abf",
  "PP.H1.abf",
  "PP.H0.abf",
  "coloc_variant_id",
  "coloc_region"
)
out_df = data.frame(matrix(vector(), 0, length(out_cols), dimnames = list(c(), out_cols)), stringsAsFactors = F)

################################################################################
if (file.size(tophits_tsv_path) == 0L) {  # empty file
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
}

################################################################################
# Load tophits
tophits_df = read.table(tophits_tsv_path, sep = '\t', header = T)
tophits_df = tophits_df[order(tophits_df$pval, decreasing=TRUE), ]  # less significant first, a priori not effect
tophits_df = tophits_df[!duplicated(tophits_df[, c('chrom', 'pos', 'nea', 'ea')]), ]
tophits_df = tophits_df[order(tophits_df$chrom, tophits_df$pos),]
if (nrow(tophits_df) == 0) {  # not tophits
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
}
tophits_df = tophits_df %>% dplyr::rename(ref = nea, alt = ea, maf = eaf)  # rename maf column
tophits_df$variant_id = paste(tophits_df$chrom, tophits_df$pos, tophits_df$ref, tophits_df$alt, sep="_")

################################################################################
# Load eqtl permuted df
permuted_df = read.table(eqtl_permuted_path, sep = "\t", header = TRUE)
permuted_df = permuted_df[permuted_df$p_perm < 0.01,]
# stop if regulatory QTL with beta-distribution permuted p value below 0.01 (bpval <0.01), 2021.Li.Mu.GenomeBiology.impactcelltype
if (nrow(permuted_df) == 0) {
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
  }

################################################################################
# Loop over tophits
coloc_variant_id = "10_60519366_C_T"
for (coloc_variant_id in unique(tophits_df$variant_id)) {
  chrom = tophits_df[tophits_df$variant_id == coloc_variant_id, "chrom"]
  pos = tophits_df[tophits_df$variant_id == coloc_variant_id, "pos"]
  start = tophits_df[tophits_df$variant_id == coloc_variant_id, "pos"] - window / 2
  if (start < 1) { start = 1 }
  end = tophits_df[tophits_df$variant_id == coloc_variant_id, "pos"] + window / 2 - 1
  if (chrom==6 & pos >= 25000000 & pos <= 35000000) { next }  # continue if MHC locus
  coloc_lead_region = paste0(chrom, ":", start, "-", end)
  print(sprintf("%s %s %s", gwas_id, eqtl_id, coloc_variant_id))

  ################################################################################
  # Load gwas summary statistics
  gwas_vcf <- gwasvcf::query_chrompos_file(chrompos = coloc_lead_region, vcffile = gwas_vcf_path)
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_vcf = gwasvcf::query_pval_vcf(vcf = gwas_vcf, pval = pval)  # select associations
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_tbl <- gwas_vcf %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()

  gwas_tbl = gwas_tbl %>% dplyr::rename(chrom = seqnames, pos = start, rsid = ID, ref = REF, alt = ALT, gwas_beta = ES, gwas_maf = AF, gwas_id=id, gwas_ss=SS, gwas_se=SE)  # rename columns
  gwas_tbl$gwas_pval = 10^-gwas_tbl$LP
  col_select=c('chrom', 'pos', 'rsid' , 'ref', 'alt', 'gwas_pval', 'gwas_beta', 'gwas_maf', 'gwas_se', 'gwas_ss')
  gwas_tbl = gwas_tbl[, col_select]
  gwas_tbl$variant_id = paste(gwas_tbl$chrom, gwas_tbl$pos, gwas_tbl$ref, gwas_tbl$alt, sep="_")
  gwas_tbl = gwas_tbl[!duplicated(gwas_tbl$variant_id),]  # keep unique RSIDs
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$gwas_beta),]  # keep non-null effect size/Z
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$gwas_se),]  # remove non-null se
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$gwas_ss),]  # keep non-null ss
  if (nrow(gwas_tbl) == 0) { next }  # continue if empty gwas

  # Load eqtl permuted and intersect gwas
  permuted_chrom_df = permuted_df[permuted_df$chrom == chrom,]
  permuted_start_df = permuted_chrom_df[permuted_chrom_df$pos >= start,]
  permuted_region_df = permuted_start_df[permuted_start_df$pos <= end,]

  ################################################################################
  # Load eqtls summary statistics
  eqtl_tbl = seqminer::tabix.read.table(tabixFile = eqtl_all_path, tabixRange = coloc_lead_region, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty
  # rename columns
  eqtl_cols = c("molecular_trait_id", "chrom", "pos", "ref", "alt", "variant_id", "ma_samples", "eqtl_maf", "eqtl_pval", "eqtl_beta", "eqtl_se", "type", "ac", "eqtl_an", "r2", "molecular_trait_object_id", "eqtl_gene_id", "median_tpm", "rsid")
  colnames(eqtl_tbl) <- eqtl_cols
  col_sel = c("molecular_trait_id", "chrom", "pos", "ref", "alt", "variant_id", "eqtl_maf", "eqtl_pval", "eqtl_beta", "eqtl_se", "eqtl_an", "eqtl_gene_id", "rsid")
  eqtl_tbl = eqtl_tbl[, col_sel]
  eqtl_tbl$variant_id = gsub("chr", "", eqtl_tbl$variant_id)  # strip newlines
  eqtl_tbl$rsid = gsub("[\r\n]", "", eqtl_tbl$rsid)  # strip newlines

  eqtl_tbl = eqtl_tbl[eqtl_tbl$eqtl_maf<1, ]  # keep MAF<1
  eqtl_tbl = eqtl_tbl[eqtl_tbl$eqtl_maf>0, ]  # keep MAF>0
  # eqtl_tbl = eqtl_tbl[!duplicated(eqtl_tbl$rsid), ]  # keep non-duplicated RSIDs
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$eqtl_beta),]  # keep non-null effect size, z, beta
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$eqtl_se),]  # remove non-null se
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty gwas

  ####################### Loop over moleculear trait id
  molecular_trait_id_lst = unique(permuted_region_df$molecular_trait_id)
  molecular_trait_id = "ILMN_2390609"
  for (molecular_trait_id in molecular_trait_id_lst) {
    print(sprintf("%s %s %s %s", gwas_id, eqtl_id, coloc_variant_id, molecular_trait_id))
    eqtl_molecular_trait_id_tbl = eqtl_tbl[eqtl_tbl$molecular_trait_id == molecular_trait_id,]

    if (nrow(eqtl_molecular_trait_id_tbl) == 0) { next }  # continue if empty gwas

    ################################################################################
    # Coloc

    # merge gwas and eqtl
    merge_df = merge(gwas_tbl,
                     eqtl_molecular_trait_id_tbl,
                     by = c("chrom", "pos", "ref", "alt", "rsid", "variant_id"))

    if (nrow(merge_df) == 0) {
      next
    }  # continue if empty merge

    # Update MAF with 1000 genomes when all NA in opengwas
    if (any(is.na(merge_df$gwas_maf))) {
      con <- dbConnect(RSQLite::SQLite(), eur_af_sqlite)
      query <- dbSendQuery(con, 'SELECT chrom, id, ref, alt, eur_af FROM eur_af where chrom=:chrom and id=:rsid and ref=:ref and alt=:alt')
      query_values = as.list(merge_df[is.na(merge_df$gwas_maf), c('chrom', 'rsid', 'ref', 'alt')])
      dbBind(query, params = query_values)
      fetch_missing_mafs_df = dbFetch(query)
      fetch_missing_mafs_df = fetch_missing_mafs_df %>% dplyr::rename(rsid = id, gwas_maf=eur_af)
      
      merge_with_maf_na_df = merge_df[is.na(merge_df$gwas_maf), c('chrom', 'pos', 'ref', 'alt', 'rsid', 'variant_id', 'gwas_pval', 'gwas_beta', 'gwas_se', 'gwas_ss', 'molecular_trait_id', 'eqtl_maf', 'eqtl_pval', 'eqtl_beta', 'eqtl_se', 'eqtl_an', 'eqtl_gene_id')]
      merge_with_maf_ok_df = merge_df[!is.na(merge_df$gwas_maf), ]
      
      merge_with_maf_new_df = merge(merge_with_maf_na_df, fetch_missing_mafs_df, by=c('chrom', 'rsid', 'ref', 'alt'))
      merge_with_maf_new_df = merge_with_maf_new_df[, colnames(merge_with_maf_ok_df)]
      
      merge_df = unique(rbind(merge_with_maf_ok_df, merge_with_maf_new_df))
    }
    if (nrow(merge_df) == 0) { next }  # continue if empty merge

    merge_df = merge_df[merge_df$gwas_maf<1, ]  # keep MAF<1
    merge_df = merge_df[merge_df$gwas_maf>0, ]  # keep MAF>0
    if (nrow(merge_df) == 0) { next }  # continue if empty merge

    # Format for coloc
    type1 = 'quant'

    gwas_coloc_lst = list(
        pvalues = merge_df$gwas_pval,
        N = merge_df$gwas_ss,
        MAF = merge_df$gwas_maf,
        beta = merge_df$gwas_beta,
        varbeta = merge_df$gwas_se ^ 2,
        type = type1,
        snp = merge_df$variant_id,
        z = merge_df$gwas_beta / merge_df$gwas_se,
        id = gwas_id
      )
    
    eqtl_coloc_lst = list(
      pvalues = merge_df$eqtl_pval,
      N = (merge_df$eqtl_an)[1] / 2,
      # Samples size is allele number (AN) dvided by 2
      MAF = merge_df$eqtl_maf,
      beta = merge_df$eqtl_beta,
      varbeta = merge_df$eqtl_se ^ 2,
      type = type1,
      snp = merge_df$variant_id,
      z = merge_df$eqtl_beta / merge_df$eqtl_se,
      id = eqtl_id
    )

    options(warn = -1)  # turn off warning
    invisible(capture.output(
      coloc_res <- coloc::coloc.abf(gwas_coloc_lst, eqtl_coloc_lst)
    ))
    options(warn = 0)  # turn on warning
    # print(coloc_res)
    ################################################################################
    # Format output
    
    coloc_res_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
    
    # coloc_res_df = coloc_res_df[coloc_df$SNP.PP.H4 > 0.2,]
    coloc_res_df = dplyr::rename(coloc_res_df, variant_id = snp)
    
    coloc_res_df = coloc_res_df %>% dplyr::mutate(nsnps = coloc_res$summary[['nsnps']])
    coloc_res_df = coloc_res_df %>% dplyr::mutate(PP.H4.abf = coloc_res$summary[['PP.H4.abf']])
    coloc_res_df = coloc_res_df %>% dplyr::mutate(PP.H3.abf = coloc_res$summary[['PP.H3.abf']])
    coloc_res_df = coloc_res_df %>% dplyr::mutate(PP.H2.abf = coloc_res$summary[['PP.H2.abf']])
    coloc_res_df = coloc_res_df %>% dplyr::mutate(PP.H1.abf = coloc_res$summary[['PP.H1.abf']])
    coloc_res_df = coloc_res_df %>% dplyr::mutate(PP.H0.abf = coloc_res$summary[['PP.H0.abf']])
    if (nrow(coloc_res_df[coloc_res_df$PP.H4.abf>=0.8, ]) > 0) {
      print( sprintf( "%s %s %s %s %s nb colocs: %d", gwas_id, eqtl_id, coloc_variant_id, coloc_lead_region, molecular_trait_id, nrow(coloc_res_df[coloc_res_df$PP.H4.abf>=0.8, ])) )
    }
    
    # print(237)
    # merge coloc results and input
    cols_select = c(
      "chrom",
      "pos",
      "rsid",
      "ref",
      "alt",
      "gwas_beta",
      "gwas_pval",
      "eqtl_beta",
      "eqtl_pval",
      'variant_id'
    )
    snp_info_df = merge_df[, cols_select]
    out_coloc_df = merge(snp_info_df, coloc_res_df, by = "variant_id")

    # add columns
    out_coloc_df$molecular_trait_id = molecular_trait_id
    out_coloc_df$eqtl_gene_id = unique(eqtl_molecular_trait_id_tbl$eqtl_gene_id)
    out_coloc_df$gwas_id = gwas_id
    out_coloc_df$eqtl_id = eqtl_id
    out_coloc_df$coloc_variant_id = coloc_variant_id
    out_coloc_df$coloc_region = coloc_lead_region
    # coloc_df$gwas_pvalue = exp(-coloc_df$gwas_pvalue)  # change -log10 pval to pval
    out_coloc_df = out_coloc_df[, out_cols]
    # print(dim(coloc_df))
    out_df = rbind(out_df, out_coloc_df)
  }
}
# print(260)
write.table(
  out_df,
  out_tsv_path,
  quote = F,
  sep = "\t",
  row.names = F,
  append = F,
  col.names = T
)
