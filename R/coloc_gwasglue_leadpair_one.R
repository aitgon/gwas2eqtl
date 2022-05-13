source("R/import_eqtl_catalogue.R")

coloc_gwasglue_leadpair_one <- function(coloc_window, eqtl_identifier, lead_egene, lead_chrom, lead_pos_hg38, eqtl_tsv_gz_path, gwas_identifier, gwas_trait_name, gwas_hg38_vcf_gz_path) {

# out_columns = c('snp', 'SNP.PP.H4', 'nsnps', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf')
# out_df <- data.frame(matrix(ncol = length(out_columns), nrow = 0))
# colnames(out_df) <- out_columns

################################################################################
#
# leadpairs
#
################################################################################

lead_region_start = as.integer(max(0, lead_pos_hg38 - coloc_window))
lead_region_end = as.integer(lead_pos_hg38 + coloc_window)
chrom_start_end_str = chrompositionstartend_str = paste0(lead_chrom, ":", lead_region_start, "-", lead_region_end)

###############################################################################
#
# eqtl
#
###############################################################################

eqtl_sumstats_df <- import_eqtl_catalogue(eqtl_sumstats_path = eqtl_tsv_gz_path,
                                          region = chrom_start_end_str,
                                          molecular_trait_id = lead_egene)

# remove duplicated snps
eqtl_sumstats_df = eqtl_sumstats_df[!duplicated(eqtl_sumstats_df$rsid), ]

if (nrow(eqtl_sumstats_df) == 0) {  # If df empty, stop and write empty file
  cat("No significant eQTLs at this p-value\n")
  return(NULL)
}


eqtl_coloc_lst = list(varbeta = eqtl_sumstats_df$se^2, 
                      N = (eqtl_sumstats_df$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                      MAF = eqtl_sumstats_df$maf, 
                      type = "quant", 
                      beta = eqtl_sumstats_df$beta,
                      snp = eqtl_sumstats_df$rsid)

################################################################################
#
# gwas
#
################################################################################

gwasr_coloc_lst <- gwasglue::gwasvcf_to_coloc(gwas_hg38_vcf_gz_path, gwas_hg38_vcf_gz_path, chrom_start_end_str)[[1]]

if (sum(gwasr_coloc_lst$snp %in% eqtl_coloc_lst$snp) == 0) {
  cat("No common SNPs between GWAS and eQTL locus\n")
  return(NULL)
}

# remove MAF==1
maf_valid = gwasr_coloc_lst$MAF<1

gwasr_coloc_lst$pvalues = gwasr_coloc_lst$pvalues[maf_valid]
gwasr_coloc_lst$N = gwasr_coloc_lst$N[maf_valid]
gwasr_coloc_lst$MAF = gwasr_coloc_lst$MAF[maf_valid]
gwasr_coloc_lst$beta = gwasr_coloc_lst$beta[maf_valid]
gwasr_coloc_lst$varbeta = gwasr_coloc_lst$varbeta[maf_valid]
gwasr_coloc_lst$snp = gwasr_coloc_lst$snp[maf_valid]
gwasr_coloc_lst$z = gwasr_coloc_lst$z[maf_valid]
gwasr_coloc_lst$chr = gwasr_coloc_lst$chr[maf_valid]
gwasr_coloc_lst$pos = gwasr_coloc_lst$pos[maf_valid]

# remove MAF==0
maf_valid = gwasr_coloc_lst$MAF>0

gwasr_coloc_lst$pvalues = gwasr_coloc_lst$pvalues[maf_valid]
gwasr_coloc_lst$N = gwasr_coloc_lst$N[maf_valid]
gwasr_coloc_lst$MAF = gwasr_coloc_lst$MAF[maf_valid]
gwasr_coloc_lst$beta = gwasr_coloc_lst$beta[maf_valid]
gwasr_coloc_lst$varbeta = gwasr_coloc_lst$varbeta[maf_valid]
gwasr_coloc_lst$snp = gwasr_coloc_lst$snp[maf_valid]
gwasr_coloc_lst$z = gwasr_coloc_lst$z[maf_valid]
gwasr_coloc_lst$chr = gwasr_coloc_lst$chr[maf_valid]
gwasr_coloc_lst$pos = gwasr_coloc_lst$pos[maf_valid]

# remove duplicated SNPs
unique_snps_bool_lst = !duplicated(gwasr_coloc_lst$snp)
gwasr_coloc_lst$pvalues = gwasr_coloc_lst$pvalues[unique_snps_bool_lst]
gwasr_coloc_lst$N = gwasr_coloc_lst$N[unique_snps_bool_lst]
gwasr_coloc_lst$MAF = gwasr_coloc_lst$MAF[unique_snps_bool_lst]
gwasr_coloc_lst$beta = gwasr_coloc_lst$beta[unique_snps_bool_lst]
gwasr_coloc_lst$varbeta = gwasr_coloc_lst$varbeta[unique_snps_bool_lst]
gwasr_coloc_lst$snp = gwasr_coloc_lst$snp[unique_snps_bool_lst]
gwasr_coloc_lst$z = gwasr_coloc_lst$z[unique_snps_bool_lst]
gwasr_coloc_lst$chr = gwasr_coloc_lst$chr[unique_snps_bool_lst]
gwasr_coloc_lst$pos = gwasr_coloc_lst$pos[unique_snps_bool_lst]

# remove z missing values
non_na_z_lst = !is.na(gwasr_coloc_lst$z)
gwasr_coloc_lst$pvalues = gwasr_coloc_lst$pvalues[non_na_z_lst]
gwasr_coloc_lst$N = gwasr_coloc_lst$N[non_na_z_lst]
gwasr_coloc_lst$MAF = gwasr_coloc_lst$MAF[non_na_z_lst]
gwasr_coloc_lst$beta = gwasr_coloc_lst$beta[non_na_z_lst]
gwasr_coloc_lst$varbeta = gwasr_coloc_lst$varbeta[non_na_z_lst]
gwasr_coloc_lst$snp = gwasr_coloc_lst$snp[non_na_z_lst]
gwasr_coloc_lst$z = gwasr_coloc_lst$z[non_na_z_lst]
gwasr_coloc_lst$chr = gwasr_coloc_lst$chr[non_na_z_lst]
gwasr_coloc_lst$pos = gwasr_coloc_lst$pos[non_na_z_lst]

################################################################################
#
# coloc.abf
#
################################################################################

coloc_res <- coloc::coloc.abf(eqtl_coloc_lst, gwasr_coloc_lst)
out_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
out_df = out_df[out_df$SNP.PP.H4>0.2, ]

out_df = out_df %>% dplyr::mutate(PP.H4.abf=coloc_res$summary[['PP.H4.abf']])
out_df = out_df %>% dplyr::mutate(nsnps=coloc_res$summary[['nsnps']])
out_df = out_df %>% dplyr::mutate(PP.H3.abf=coloc_res$summary[['PP.H3.abf']])
out_df = out_df %>% dplyr::mutate(PP.H2.abf=coloc_res$summary[['PP.H2.abf']])
out_df = out_df %>% dplyr::mutate(PP.H1.abf=coloc_res$summary[['PP.H1.abf']])
out_df = out_df %>% dplyr::mutate(PP.H0.abf=coloc_res$summary[['PP.H0.abf']])

out_df = out_df %>% dplyr::mutate(lead_chrom = lead_chrom)
out_df = out_df %>% dplyr::mutate(lead_pos_hg38 = lead_pos_hg38)
out_df = out_df %>% dplyr::mutate(lead_gene = lead_egene)
out_df = out_df %>% dplyr::mutate(eqtl_identifier = eqtl_identifier)
out_df = out_df %>% dplyr::mutate(gwas_identifier = gwas_identifier)
out_df = out_df %>% dplyr::mutate(gwas_trait_name = gwas_trait_name)

###############################################################################
#
# Add eQTL information to colocalized SNPs
# Format output table
#
###############################################################################

# eqtl_sumstats_df = dplyr::rename(eqtl_sumstats_df, rsid=rsid)
out_df = dplyr::rename(out_df, rsid=snp)


  # Merge
out_df = merge(out_df, eqtl_sumstats_df[, c("rsid", "position", "molecular_trait_object_id", "beta", "pvalue", "ref", "alt")], by="rsid")

  # Rename columns

out_df = dplyr::rename(out_df, pos=position)
out_df = dplyr::rename(out_df, ref=ref)
out_df = dplyr::rename(out_df, alt=alt)
out_df = dplyr::rename(out_df, egene=molecular_trait_object_id)
out_df = dplyr::rename(out_df, eqtl_beta=beta)
out_df = dplyr::rename(out_df, eqtl_pvalue=pvalue)
out_df = dplyr::rename(out_df, pp_h4=SNP.PP.H4)

out_df = dplyr::rename(out_df, chrom=lead_chrom)
out_df = dplyr::rename(out_df, leadeqtl_pos=lead_pos_hg38)
out_df = dplyr::rename(out_df, leadeqtl_egene=lead_gene)

################################################################################
#
# Add GWAS information to colocalized SNPs
# Format output table
#
################################################################################


gwas_df = bind_rows(gwasr_coloc_lst)
gwas_df = dplyr::rename(gwas_df, chrom=chr)
gwas_df = dplyr::rename(gwas_df, rsid=snp)
gwas_df = dplyr::rename(gwas_df, pos=pos)

out_df = merge(out_df, gwas_df[, c("chrom", "rsid", "pos", "beta", "pvalues")], by=c("chrom", "rsid", "pos"))

out_df = dplyr::rename(out_df, gwas_beta=beta)
out_df = dplyr::rename(out_df, gwas_pvalue=pvalues)

out_df = out_df %>% select(chrom, pos, rsid, ref, alt, egene, eqtl_beta, eqtl_pvalue, eqtl_identifier, gwas_beta, gwas_pvalue, gwas_identifier, gwas_trait_name, pp_h4, PP.H4.abf, nsnps, PP.H3.abf, PP.H2.abf, PP.H1.abf, PP.H0.abf, leadeqtl_pos, leadeqtl_egene)
#print(out_df)
return(out_df)
}
