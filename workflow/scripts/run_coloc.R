suppressWarnings(suppressPackageStartupMessages({
library(VariantAnnotation)
library(dplyr)
}))

options(warn=-1)  # turn off warning

# window=500000
# pcutoff=5e-8
# gwas_vcf_path = "/home/gonzalez/Software/process/hg38/gwas.mrcieu.ac.uk/files/ieu-a-1162/ieu-a-1162.vcf.bgz"
# eqtl_leads_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.leadpair.tsv"
# eqtl_all_path = "/home/gonzalez/Software/public/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.all.tsv.gz"
# out_tsv_path = "coloc.tsv"

args = commandArgs(trailingOnly=TRUE)
window = as.numeric(args[1])
pcutoff = as.numeric(args[2])
gwas_vcf_path = args[3]
eqtl_leads_path = args[4]
eqtl_all_path = args[5]
out_tsv_path = args[6]

eqtl_identifier = gsub(".all.tsv.gz", "", strsplit(eqtl_all_path, split="/", fixed=T)[[1]][length(strsplit(eqtl_all_path, split="/", fixed=T)[[1]])])
gwas_identifier = strsplit(gwas_vcf_path, split="/", fixed=T)[[1]][length(strsplit(gwas_vcf_path, split="/", fixed=T)[[1]])-1]
dir.create(dirname(out_tsv_path), showWarnings = FALSE, recursive=TRUE)

# no bcftools error message and exit
if (length(which(!is.na(system('which bcftools', intern=T))))==0) {
  stop("Error. Please, install bcftools")
} else {
  bcftools_path = system("which bcftools", intern=T)
  gwasvcf::set_bcftools(bcftools_path)
}

################################################################################
# Init output df
mat = matrix(ncol = 0, nrow = 0)
out_columns = c("chrom ", "pos ", "rsid ", "ref ", "alt ", "egene ", "eqtl_beta ", "eqtl_pvalue ", "eqtl_identifier ", "gwas_beta ", "gwas_pvalue ", "gwas_identifier ", "pp_h4 ", "PP.H4.abf ", "coloc_window ", "nsnps ", "PP.H3.abf ", "PP.H2.abf ", "PP.H1.abf ", "PP.H0.abf")
out_df = data.frame(matrix(vector(), 0, length(out_columns), dimnames=list(c(), out_columns)),  stringsAsFactors=F)

################################################################################
# Load lead eqtl pairs
eqtl_leads_df = read.table(eqtl_leads_path, sep='\t', header=T)
eqtl_leads_df = eqtl_leads_df[eqtl_leads_df$chromosome %in% as.character(1:22), ]
eqtl_leads_df$start = eqtl_leads_df$position - (window/2)
eqtl_leads_df$end = eqtl_leads_df$position + (window/2)
eqtl_leads_df[eqtl_leads_df$start<0, "start"] <-0  # negative values with 0
# region_lst = as.character(paste0(eqtl_leads_df$chrom, ':', eqtl_leads_df$start, '-', eqtl_leads_df$end))

################################################################################
# Loop over regions
# region_i = region_lst[1]
# region_i = "1:120860278-121860278"
# i = 86
for (i in 1:nrow(eqtl_leads_df)) {
  
  chrom = eqtl_leads_df[i, "chromosome"]
  start = eqtl_leads_df[i, "start"]
  end = eqtl_leads_df[i, "end"]
  molecular_trait_id_i = eqtl_leads_df[i, "molecular_trait_id"]
  region_i = paste0(chrom, ":", start, "-", end)
  print(region_i)

  ################################################################################
  # Load and filter gwas summary statistics
  gwas_vcf <- gwasvcf::query_chrompos_file(chrompos=region_i, vcffile=gwas_vcf_path)
  if(length(gwas_vcf) == 0) {next}  # continue if empty gwas
  gwas_vcf <- gwasvcf::query_pval_vcf(vcf=gwas_vcf, pval=pcutoff)  # select associations
  if(length(gwas_vcf) == 0) {next}  # continue if empty gwas
  gwas_tbl <- gwas_vcf %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  
  gwas_tbl = gwas_tbl[gwas_tbl$AF<1, ]  # keep MAF<1
  gwas_tbl = gwas_tbl[gwas_tbl$AF>0, ]  # keep MAF>0
  gwas_tbl = gwas_tbl[!duplicated(gwas_tbl$ID), ]  # keep unique RSIDs
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$ES), ]  # keep non-null effect size/Z
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SE), ]  # remove non-null se
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SS), ]  # keep non-null ss
  if(nrow(gwas_tbl) == 0) {next}  # continue if empty gwas
  
  ################################################################################
  # Load eqtls summary statistics
  eqtl_tbl = seqminer::tabix.read.table(tabixFile = eqtl_all_path, tabixRange = region_i, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  if(nrow(eqtl_tbl) == 0) {next}  # continue if empty gwas
  
  eqtl_cols = c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
  colnames(eqtl_tbl) <- eqtl_cols
  eqtl_tbl$rsid = gsub("\r", "", eqtl_tbl$rsid) # strip \r from rsid
  
  eqtl_tbl = eqtl_tbl[eqtl_tbl$maf<1, ]  # keep MAF<1
  eqtl_tbl = eqtl_tbl[eqtl_tbl$maf>0, ]  # keep MAF>0
  # eqtl_tbl = eqtl_tbl[!duplicated(eqtl_tbl$rsid), ]  # keep non-duplicated RSIDs
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$beta), ]  # keep non-null effect size, z, beta
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$se), ]  # remove non-null se
  if(nrow(eqtl_tbl) == 0) {next}  # continue if empty gwas
  
  ################################################################################
  # Get eqtl_egene 
    
    # print(molecular_trait_id_i)
    
    egene_ensg = as.character(eqtl_tbl[head(which(eqtl_tbl$molecular_trait_id==molecular_trait_id_i), 1), "molecular_trait_object_id"])
    # print(egene_ensg)

    eqtl_egene_tbl = eqtl_tbl[eqtl_tbl$molecular_trait_id==molecular_trait_id_i, ]
    if(nrow(eqtl_egene_tbl) == 0) {next}  # continue if empty gwas
    
    ################################################################################
    # Coloc
    
    # Keep common SNPs
    merge_df = merge(eqtl_egene_tbl, gwas_tbl, by.x=c("rsid", "ref", "alt", "chromosome", "position"), by.y=c("ID", "REF", "ALT", "seqnames", "start"))
    if(nrow(merge_df) == 0) {next}  # continue if empty gwas
    # gwas_tbl = gwas_tbl[gwas_tbl$ID %in% rsid_intersection_lst, ]  # keep common rsids
    # eqtl_egene_tbl = eqtl_egene_tbl[eqtl_egene_tbl$rsid %in% rsid_intersection_lst, ]  # keep common rsids
    
    # Format for coloc
    type1='quant'
    eqtl_coloc_lst = list(pvalues = merge_df$pvalue,
                          N = (merge_df$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                          MAF = merge_df$maf, 
                          beta = merge_df$beta,
                          varbeta = merge_df$se^2, 
                          type = type1, 
                          snp = merge_df$rsid)
    
    gwas_coloc_lst = merge_df %>% {list(pvalues = 10^-.$LP, 
                                        N = .$SS, 
                                        MAF = .$AF, 
                                        beta = .$ES,
                                        varbeta = .$SE^2,
                                        type = type1, 
                                        snp = .$rsid, 
                                        z = .$ES / .$SE, 
                                        id = VariantAnnotation::samples(VariantAnnotation::header(gwas_vcf))[1])}
    
    options(warn=-1)  # turn off warning
    invisible(capture.output(coloc_res <- coloc::coloc.abf(eqtl_coloc_lst, gwas_coloc_lst)))
    options(warn=0)  # turn on warning
    
    ################################################################################
    # Format output
    
    coloc_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
    coloc_df = coloc_df[coloc_df$SNP.PP.H4>0.2, ]
    
    coloc_df = coloc_df %>% dplyr::mutate(PP.H4.abf=coloc_res$summary[['PP.H4.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(nsnps=coloc_res$summary[['nsnps']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H3.abf=coloc_res$summary[['PP.H3.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H2.abf=coloc_res$summary[['PP.H2.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H1.abf=coloc_res$summary[['PP.H1.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H0.abf=coloc_res$summary[['PP.H0.abf']])
    
    # merge coloc results and input
    coloc_df = merge(coloc_df, merge_df, by.x="snp", by.y="rsid")
    
    # rename columns
    coloc_df = dplyr::rename(coloc_df, chrom=chromosome, pos=position, rsid=snp, egene=gene_id, eqtl_beta=beta, eqtl_pvalue=pvalue, gwas_beta=ES, gwas_pvalue=LP, pp_h4=SNP.PP.H4)
    # add columns
    coloc_df = coloc_df %>% dplyr::mutate(gwas_identifier=gwas_identifier)
    coloc_df = coloc_df %>% dplyr::mutate(eqtl_identifier=eqtl_identifier)
    coloc_df = coloc_df %>% dplyr::mutate(coloc_window=region_i)
    coloc_df$gwas_pvalue = exp(-coloc_df$gwas_pvalue)  # change -log10 pval to pval
    coloc_df = coloc_df[, c("chrom", "pos", "rsid", "ref", "alt", "egene", 
                            "eqtl_beta", "eqtl_pvalue", "eqtl_identifier", 
                            "gwas_beta", "gwas_pvalue", "gwas_identifier", 
                            "pp_h4", "PP.H4.abf", "coloc_window", "nsnps", 
                            "PP.H3.abf", "PP.H2.abf", "PP.H1.abf", "PP.H0.abf")]

    out_df = rbind(out_df, coloc_df)
}

write.table(out_df, out_tsv_path, quote=F, sep="\t", row.names = F, append=F, col.names = T)
