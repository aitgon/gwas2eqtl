suppressWarnings(suppressPackageStartupMessages({
  library(DBI)
  library(VariantAnnotation)
  library(dplyr)
}))

# PARAMS
# window = 1000000
# pval = 5e-08
# r2 = 0.1
# kb = 1000
# tophits_tsv_path = sprintf("/home/gonzalez/Repositories/eqtl2gwas/out/maf/tophits/%s/pval_%s/r2_%.1f/kb_%d/hg38.tsv", gwas_id, as.character(pval), r2, kb)
# gwas_vcf_path = "/home/gonzalez/Software/process/hg38/gwas.mrcieu.ac.uk/files/ebi-a-GCST002318/ebi-a-GCST002318.vcf.bgz"
#eqtl_permuted_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.leadpair.tsv"
# eqtl_all_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_naive.all.tsv.gz"
# out_tsv_path = "coloc.tsv"
# eur_af_sqlite = "/home/gonzalez/Repositories/eqtl2gwas/out/maf/eur_af.sqlite"
# eqtl_permuted_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_naive.permuted.tsv.gz"
# END PARAMS

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop("8 arguments must be supplied", call.=FALSE)
}
window = as.numeric(args[1])
pval = as.numeric(args[2])
tophits_tsv_path = args[3]
gwas_vcf_path = args[4]
eqtl_permuted_path = args[5]
eqtl_all_path = args[6]
eur_af_sqlite = args[7]
out_tsv_path = args[8]

eqtl_id = gsub(".all.tsv.gz", "", strsplit(eqtl_all_path, split = "/", fixed = T)[[1]][length(strsplit(eqtl_all_path, split = "/", fixed = T)[[1]])])
gwas_id = strsplit(gwas_vcf_path, split = "/", fixed = T)[[1]][length(strsplit(gwas_vcf_path, split = "/", fixed = T)[[1]]) - 1]

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
out_cols = c("chrom", "pos", "rsid", "ref", "alt", "egene",
                            "gwas_beta", "gwas_pval", "gwas_id",
                            "eqtl_beta", "eqtl_pval", "eqtl_id",
                            "PP.H4.abf", "SNP.PP.H4", "nsnps",
                            "PP.H3.abf", "PP.H2.abf", "PP.H1.abf", "PP.H0.abf",
                            "coloc_lead_pos", "coloc_lead_rsid", "coloc_region")
out_df = data.frame(matrix(vector(), 0, length(out_cols), dimnames = list(c(), out_cols)), stringsAsFactors = F)

################################################################################
if (file.size(tophits_tsv_path) == 0L) {  # empty file
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
}

# Load tophits
tophits_df = read.table(tophits_tsv_path, sep = '\t', header = T)
tophits_df = tophits_df[order(tophits_df$p), ]  # most significant first 
tophits_df = tophits_df[!duplicated(tophits_df[, c('chr', 'position', 'rsid', 'nea', 'ea')]), ]

if (nrow(tophits_df) == 0) {  # not tophits
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
}

tophits_df = tophits_df %>% dplyr::rename(chrom = chr, pos = position, ref = nea, alt = ea, maf = eaf)  # rename maf column

################################################################################
# Loop over tophits
# rsid = "rs28411352"
# rsid = "rs6930468"
# rsid = "rs3130663"
# rsid_lst = c("rs9357094", "rs112702727", "rs17875360", "rs4711222", "rs3130663")
for (rsid in unique(tophits_df$rsid)) {
  # print(rsid)
  chrom = tophits_df[tophits_df$rsid == rsid, "chrom"]
  pos = tophits_df[tophits_df$rsid == rsid, "pos"]
  start = tophits_df[tophits_df$rsid == rsid, "pos"] - window / 2
  if (start < 1) { start = 1 }
  end = tophits_df[tophits_df$rsid == rsid, "pos"] + window / 2
  if (chrom==6 & pos >= 25000000 & pos <= 35000000) { next }  # continue if MHC locus
  region_i = paste0(chrom, ":", start, "-", end)
  # print(paste(region_i, rsid))

  ################################################################################
  # Load gwas summary statistics
  gwas_vcf <- gwasvcf::query_chrompos_file(chrompos = region_i, vcffile = gwas_vcf_path)
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_vcf = gwasvcf::query_pval_vcf(vcf = gwas_vcf, pval = pval)  # select associations
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_tbl <- gwas_vcf %>%
    gwasvcf::vcf_to_granges() %>%
    dplyr::as_tibble()

  gwas_tbl = gwas_tbl %>% dplyr::rename(chrom = seqnames, pos = start, rsid = ID, ref = REF, alt = ALT, maf = AF)  # rename columns
  gwas_tbl$gwas_pval = 10^-gwas_tbl$LP
  gwas_tbl = gwas_tbl %>% dplyr::rename(gwas_beta = ES, gwas_maf = maf)  # rename
  # gwas_tbl = gwas_tbl[gwas_tbl$AF<1, ]  # keep MAF<1
  # gwas_tbl = gwas_tbl[gwas_tbl$AF>0, ]  # keep MAF>0
  gwas_tbl = gwas_tbl[!duplicated(gwas_tbl$rsid),]  # keep unique RSIDs
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$gwas_beta),]  # keep non-null effect size/Z
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SE),]  # remove non-null se
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SS),]  # keep non-null ss
  if (nrow(gwas_tbl) == 0) { next }  # continue if empty gwas

  # Load eqtl permuted and intersect gwas
  permuted_df = read.table(eqtl_permuted_path, sep = "\t", header = TRUE)
  permuted_df = permuted_df %>% dplyr::rename(chrom = chromosome, pos = position, egene = molecular_trait_id)
  permuted_df = permuted_df[permuted_df$p_perm < 0.01,]
  gwas_eqtl_permuted_df = merge(gwas_tbl, permuted_df, by = c("chrom", "pos"))
  # continue if regulatory QTL with beta-distribution permuted p value below 0.01 (bpval <0.01), 2021.Li.Mu.GenomeBiology.impactcelltype
  if (nrow(gwas_eqtl_permuted_df) == 0) { next }

  ################################################################################
  # Load eqtls summary statistics
  eqtl_tbl = seqminer::tabix.read.table(tabixFile = eqtl_all_path, tabixRange = region_i, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty
  # rename columns
  eqtl_cols = c("egene", "chrom", "pos", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
  colnames(eqtl_tbl) <- eqtl_cols
  # strip newlines
  eqtl_tbl$rsid = gsub("[\r\n]", "", eqtl_tbl$rsid)

  eqtl_tbl = eqtl_tbl[eqtl_tbl$maf<1, ]  # keep MAF<1
  eqtl_tbl = eqtl_tbl[eqtl_tbl$maf>0, ]  # keep MAF>0
  # eqtl_tbl = eqtl_tbl[!duplicated(eqtl_tbl$rsid), ]  # keep non-duplicated RSIDs
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$beta),]  # keep non-null effect size, z, beta
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$se),]  # remove non-null se
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty gwas

  ####################### Loop over moleculear trait id
  egene_lst = unique(gwas_eqtl_permuted_df$egene)
  # egene="ENSG00000204084"
  # egene = "ENSG00000204528"
  for (egene in egene_lst) {
    # print(egene)

    eqtl_egene_tbl = eqtl_tbl[eqtl_tbl$egene == egene, ]
    if (nrow(eqtl_egene_tbl) == 0) { next }  # continue if empty gwas

    ################################################################################
    # Coloc

    eqtl_egene_tbl = eqtl_egene_tbl %>% dplyr::rename(eqtl_beta = beta, eqtl_pval = pvalue, eqtl_maf = maf)  # rename

    # Keep common SNPs
    merge_df = merge(gwas_tbl, eqtl_egene_tbl, by = c("chrom", "pos", "rsid", "ref", "alt"))

    if (nrow(merge_df) == 0) { next }  # continue if empty merge

    # Update MAF with 1000 genomes when all NA in opengwas
    if (all(is.na(merge_df$"gwas_maf"))) {
      # gwas_maf_na_df = merge_df[is.na(merge_df$AF), c("chromosome", "rsid", "ref", "alt", "AF")]
      con <- dbConnect(RSQLite::SQLite(), eur_af_sqlite)
      # dbGetQuery(con, 'SELECT * FROM eur_af where chrom=1 and id="rs12137845" and ref="T" and alt="C"')
      # query <- dbSendQuery(con, 'SELECT * FROM iris WHERE "Sepal.Length" < :x')
      query <- dbSendQuery(con, 'SELECT chrom, id, ref, alt, eur_af FROM eur_af where chrom=:chromosome and id=:rsid and ref=:ref and alt=:alt')
      # query <- dbSendQuery(con, 'SELECT * FROM eur_af where id=:x')
      dbBind(query, params = list(chromosome = merge_df$chrom,
                               rsid = merge_df$rsid,
                               ref = merge_df$ref,
                               alt = merge_df$alt
      ))
      gwas_maf_df = dbFetch(query)
      dbDisconnect(con)
      # gwas_maf_df = subset(gwas_maf_df, select = -c(pos) )  # drop pos hg19
      gwas_maf_df = gwas_maf_df %>% dplyr::rename(chrom = chrom,
                                                  rsid = id,
                                                  gwas_maf = eur_af,
      )  # rename maf column
      merge_df = subset(merge_df, select = -c(gwas_maf))  # drop old gwas MAF
      merge_df = merge(merge_df, gwas_maf_df, on = c("chromosome", "rsid", "ref", "alt"))
    }
    if (nrow(merge_df) == 0) { next }  # continue if empty merge

    merge_df = merge_df[merge_df$gwas_maf<1, ]  # keep MAF<1
    merge_df = merge_df[merge_df$gwas_maf>0, ]  # keep MAF>0
    if (nrow(merge_df) == 0) { next }  # continue if empty merge

    # Format for coloc
    type1 = 'quant'

    gwas_coloc_lst = merge_df %>% { list(pvalues = .$gwas_pval,
                                         N = .$SS,
                                         MAF = .$gwas_maf,
                                         beta = .$gwas_beta,
                                         varbeta = .$SE^2,
                                         type = type1,
                                         snp = .$rsid,
                                         z = .$ES / .$SE,
                                         id = VariantAnnotation::samples(VariantAnnotation::header(gwas_vcf))[1]) }

    eqtl_coloc_lst = list(pvalues = merge_df$eqtl_pval,
                          N = (merge_df$an)[1] / 2, # Samples size is allele number (AN) dvided by 2
                          MAF = merge_df$eqtl_maf,
                          beta = merge_df$eqtl_beta,
                          varbeta = merge_df$se^2,
                          type = type1,
                          snp = merge_df$rsid)

    options(warn = -1)  # turn off warning
    invisible(capture.output(coloc_res <- coloc::coloc.abf(gwas_coloc_lst, eqtl_coloc_lst)))
    options(warn = 0)  # turn on warning

    ################################################################################
    # Format output

    coloc_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
    # coloc_df = coloc_df[coloc_df$SNP.PP.H4 > 0.2,]
    coloc_df = dplyr::rename(coloc_df, rsid = snp)


    coloc_df = coloc_df %>% dplyr::mutate(nsnps = coloc_res$summary[['nsnps']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H4.abf = coloc_res$summary[['PP.H4.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H3.abf = coloc_res$summary[['PP.H3.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H2.abf = coloc_res$summary[['PP.H2.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H1.abf = coloc_res$summary[['PP.H1.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H0.abf = coloc_res$summary[['PP.H0.abf']])

    # merge coloc results and input
    merge_cols = c("chrom", "pos", "rsid", "ref", "alt", "egene", "gwas_beta", "gwas_pval", "eqtl_beta", "eqtl_pval")
    snp_info_df = merge_df[, merge_cols]
    # coloc_cols = c("chrom", "pos", "rsid", "ref", "alt", "egene", "SNP.PP.H4", 'PP.H4.abf', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', "nsnps")
    # coloc_df = coloc_df[, coloc_cols]
    coloc_df = merge(snp_info_df, coloc_df, by = "rsid")

    # rename columns
    # coloc_df = dplyr::rename(coloc_df, egene=gene_id, eqtl_beta=beta, eqtl_pvalue=pvalue, gwas_beta=ES, gwas_pvalue=LP)
    # add columns
    coloc_df$gwas_id = gwas_id
    coloc_df$eqtl_id = eqtl_id
    coloc_df$coloc_lead_pos = pos
    coloc_df$coloc_lead_rsid = rsid
    coloc_df$coloc_region = region_i
    # coloc_df$gwas_pvalue = exp(-coloc_df$gwas_pvalue)  # change -log10 pval to pval
    coloc_df = coloc_df[, out_cols]
# print(dim(coloc_df))
    out_df = rbind(out_df, coloc_df)
  }
}

write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
