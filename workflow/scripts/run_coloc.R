suppressWarnings(suppressPackageStartupMessages({
  library(coloc)
  library(DBI)
  library(dplyr)
  library(gwasvcf)
  library(VariantAnnotation)
}))

# PARAMS
gwas_id = "ebi-a-GCST002318"
eqtl_id = "Schmiedel_2018_ge_CD8_T-cell_naive"
window = 1000000
pval = 5e-08
tophits_tsv_path = "out/gwas418/tophits/ebi-a-GCST002318/pval_5e-08/r2_0.1/kb_1000/hg38.tsv"
gwas_vcf_path = "/home/gonzalez/Software/process/hg38/gwas.mrcieu.ac.uk/files/ebi-a-GCST002318/ebi-a-GCST002318.vcf.bgz"
eqtl_permuted_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_naive.permuted.tsv.gz"
eqtl_all_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_naive.all.tsv.gz"
eur_af_sqlite = "out/eur_af.sqlite"
out_tsv_path = "out/gwas418/coloc/ebi-a-GCST002318/pval_5e-08/r2_0.1/kb_1000/window_1000000/Schmiedel_2018_ge_CD8_T-cell_naive.tsv"
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

# eqtl_id = gsub(".all.tsv.gz", "", strsplit(eqtl_all_path, split = "/", fixed = T)[[1]][length(strsplit(eqtl_all_path, split = "/", fixed = T)[[1]])])
# gwas_id = strsplit(gwas_vcf_path, split = "/", fixed = T)[[1]][length(strsplit(gwas_vcf_path, split = "/", fixed = T)[[1]]) - 1]

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
  "gene_id",
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
  "coloc_lead_pos",
  "coloc_lead_variant",
  "coloc_region"
)
out_df = data.frame(matrix(vector(), 0, length(out_cols), dimnames = list(c(), out_cols)), stringsAsFactors = F)

################################################################################
if (file.size(tophits_tsv_path) == 0L) {  # empty file
  write.table(out_df, out_tsv_path, quote = F, sep = "\t", row.names = F, append = F, col.names = T)
  quit(save = "no", status = 0, runLast = FALSE)
}

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

################################################################################
# Loop over tophits
# coloc_lead_rsid = "rs28411352"
# rsid = "rs6930468"
rsid = "rs11203201"
# rsid_lst = c("rs9357094", "rs112702727", "rs17875360", "rs4711222", "rs3130663")
for (coloc_lead_rsid in unique(tophits_df$rsid)) {
  chrom = tophits_df[tophits_df$rsid == coloc_lead_rsid, "chrom"]
  pos = tophits_df[tophits_df$rsid == coloc_lead_rsid, "pos"]
  start = tophits_df[tophits_df$rsid == coloc_lead_rsid, "pos"] - window / 2
  if (start < 1) { start = 1 }
  end = tophits_df[tophits_df$rsid == coloc_lead_rsid, "pos"] + window / 2 - 1
  if (chrom==6 & pos >= 25000000 & pos <= 35000000) { next }  # continue if MHC locus
  coloc_lead_region = paste0(chrom, ":", start, "-", end)
  print(sprintf("%s %s %s", gwas_id, eqtl_id, coloc_lead_region))

  ################################################################################
  # Load gwas summary statistics
  gwas_vcf <- gwasvcf::query_chrompos_file(chrompos = coloc_lead_region, vcffile = gwas_vcf_path)
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_vcf = gwasvcf::query_pval_vcf(vcf = gwas_vcf, pval = pval)  # select associations
  if (length(gwas_vcf) == 0) { next }  # continue if empty gwas
  gwas_tbl <- gwas_vcf %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()

  gwas_tbl = gwas_tbl %>% dplyr::rename(chrom = seqnames, pos = start, rsid = ID, ref = REF, alt = ALT, gwas_maf = AF, gwas_beta = ES)  # rename columns
  gwas_tbl$gwas_pval = 10^-gwas_tbl$LP
  gwas_tbl = gwas_tbl[!duplicated(gwas_tbl$rsid),]  # keep unique RSIDs
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$gwas_beta),]  # keep non-null effect size/Z
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SE),]  # remove non-null se
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SS),]  # keep non-null ss
  if (nrow(gwas_tbl) == 0) { next }  # continue if empty gwas

  # Load eqtl permuted and intersect gwas
  permuted_df = read.table(eqtl_permuted_path, sep = "\t", header = TRUE)
  # permuted_df = permuted_df %>% dplyr::rename(chrom = chromosome, pos = position, egene = molecular_trait_object_id)
  # permuted_df = permuted_df[permuted_df$p_perm < 0.01,]
  # gwas_eqtl_permuted_df = merge(gwas_tbl, permuted_df, by = c("chrom", "pos"))
  gwas_eqtl_permuted_1_df = dplyr::filter(permuted_df, chrom == chrom & pos >= start & pos <= end)
  gwas_eqtl_permuted_df = gwas_eqtl_permuted_1_df[gwas_eqtl_permuted_1_df$chrom == chrom,]
  # continue if regulatory QTL with beta-distribution permuted p value below 0.01 (bpval <0.01), 2021.Li.Mu.GenomeBiology.impactcelltype
  if (nrow(gwas_eqtl_permuted_df) == 0) { next }

  ################################################################################
  # Load eqtls summary statistics
  eqtl_tbl = seqminer::tabix.read.table(tabixFile = eqtl_all_path, tabixRange = coloc_lead_region, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty
  # rename columns
  eqtl_cols = c("molecular_trait_id", "chrom", "pos", "ref", "alt", "variant", "ma_samples", "eqtl_maf", "eqtl_pval", "eqtl_beta", "se", "type", "ac", "an", "r2", "egene", "gene_id", "median_tpm", "rsid")
  colnames(eqtl_tbl) <- eqtl_cols
  # strip newlines
  eqtl_tbl$rsid = gsub("[\r\n]", "", eqtl_tbl$rsid)

  eqtl_tbl = eqtl_tbl[eqtl_tbl$eqtl_maf<1, ]  # keep MAF<1
  eqtl_tbl = eqtl_tbl[eqtl_tbl$eqtl_maf>0, ]  # keep MAF>0
  # eqtl_tbl = eqtl_tbl[!duplicated(eqtl_tbl$rsid), ]  # keep non-duplicated RSIDs
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$eqtl_beta),]  # keep non-null effect size, z, beta
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$se),]  # remove non-null se
  if (nrow(eqtl_tbl) == 0) { next }  # continue if empty gwas

  ####################### Loop over moleculear trait id
  # egene_lst = unique(gwas_eqtl_permuted_df$egene)
  molecular_trait_id_lst = unique(gwas_eqtl_permuted_df$molecular_trait_id)
  # egene="ENSG00000204084"
  molecular_trait_id = "ENSG00000160185"
  for (molecular_trait_id in molecular_trait_id_lst) {
    # print(
    #   sprintf(
    #     "%s %s %s %s %s %s",
    #     gwas_id,
    #     eqtl_id,
    #     chrom,
    #     pos,
    #     coloc_lead_rsid,
    #     molecular_trait_id
    #   )
    # )

    eqtl_molecular_trait_id_tbl = eqtl_tbl[eqtl_tbl$molecular_trait_id == molecular_trait_id,]
    
    # remove duplicate snps
    # eqtl_egene_tbl = eqtl_egene_tbl[order(eqtl_egene_tbl$pvalue, decreasing=TRUE), ]  # less significant first, a priori not effect
    # eqtl_egene_tbl = eqtl_egene_tbl[!duplicated(eqtl_egene_tbl[, c('chrom', 'pos', 'rsid', 'ref', 'alt', 'egene')]), ]
    
    if (nrow(eqtl_molecular_trait_id_tbl) == 0) { next }  # continue if empty gwas

    ################################################################################
    # Coloc

    # eqtl_egene_tbl = eqtl_egene_tbl %>% dplyr::rename(eqtl_beta = beta, eqtl_pval = pvalue, eqtl_maf = maf)  # rename

    # merge gwas and eqtl
    merge_df = merge(gwas_tbl,
                     eqtl_molecular_trait_id_tbl,
                     by = c("chrom", "pos", "ref", "alt", "rsid"))
    # create variant_id based on chrom, pos, ref and alt
    merge_df$variant_id = paste(merge_df$chrom,
                                merge_df$pos,
                                merge_df$ref,
                                merge_df$alt,
                                sep = "_")
    coloc_cols = c(
      "chrom",
      "pos",
      "rsid",
      "ref",
      "alt",
      "gwas_pval",
      "gwas_beta",
      "gwas_maf",
      "SS",
      "SE",
      "eqtl_pval",
      "eqtl_beta",
      "eqtl_maf",
      "an",
      "se",
      "variant_id"
    )
    merge_df = merge_df[, coloc_cols]
    merge_df = na.omit(merge_df)
    if (nrow(merge_df) == 0) {
      next
    }  # continue if empty merge

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

    gwas_coloc_lst = merge_df %>% {
      list(
        pvalues = .$gwas_pval,
        N = .$SS,
        MAF = .$gwas_maf,
        beta = .$gwas_beta,
        varbeta = .$SE ^ 2,
        type = type1,
        snp = .$rsid,
        z = .$gwas_beta / .$SE,
        id = gwas_id
      )
    }
    
    eqtl_coloc_lst = list(
      pvalues = merge_df$eqtl_pval,
      N = (merge_df$an)[1] / 2,
      # Samples size is allele number (AN) dvided by 2
      MAF = merge_df$eqtl_maf,
      beta = merge_df$eqtl_beta,
      varbeta = merge_df$se ^ 2,
      type = type1,
      snp = merge_df$rsid,
      id = eqtl_id
    )

    options(warn = -1)  # turn off warning
    invisible(capture.output(
      coloc_res <- coloc::coloc.abf(gwas_coloc_lst, eqtl_coloc_lst)
    ))
    options(warn = 0)  # turn on warning
    print(coloc_res)
    ################################################################################
    # Format output
    
    coloc_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
    # coloc_df = coloc_df[coloc_df$SNP.PP.H4 > 0.2,]
    coloc_df = dplyr::rename(coloc_df, variant_id = snp)
    
    coloc_df = coloc_df %>% dplyr::mutate(nsnps = coloc_res$summary[['nsnps']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H4.abf = coloc_res$summary[['PP.H4.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H3.abf = coloc_res$summary[['PP.H3.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H2.abf = coloc_res$summary[['PP.H2.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H1.abf = coloc_res$summary[['PP.H1.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H0.abf = coloc_res$summary[['PP.H0.abf']])
    # print(237)
    # merge coloc results and input
    merge_cols = c(
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
    snp_info_df = merge_df[, merge_cols]
    # coloc_cols = c("chrom", "pos", "rsid", "ref", "alt", "egene", "SNP.PP.H4", 'PP.H4.abf', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', "nsnps")
    # coloc_df = coloc_df[, coloc_cols]
    coloc_df = merge(snp_info_df, coloc_df, by = "variant_id")
    # print(244)
    # rename columns
    # coloc_df = dplyr::rename(coloc_df, egene=gene_id, eqtl_beta=beta, eqtl_pvalue=pvalue, gwas_beta=ES, gwas_pvalue=LP)
    # add columns
    coloc_df$molecular_trait_id = molecular_trait_id
    coloc_df$gene_id = unique(eqtl_molecular_trait_id_tbl$gene_id)
    coloc_df$gwas_id = gwas_id
    coloc_df$eqtl_id = eqtl_id
    coloc_df$coloc_lead_pos = pos
    coloc_df$coloc_lead_variant = coloc_lead_variant_id
    coloc_df$coloc_region = coloc_lead_region
    # coloc_df$gwas_pvalue = exp(-coloc_df$gwas_pvalue)  # change -log10 pval to pval
    coloc_df = coloc_df[, out_cols]
    # print(dim(coloc_df))
    out_df = rbind(out_df, coloc_df)
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
