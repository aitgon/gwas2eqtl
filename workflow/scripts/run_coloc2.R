suppressWarnings(suppressPackageStartupMessages({
library(VariantAnnotation)
library(dplyr)
}))

# window=500000
# pcutoff=5e-8
# gwas_vcf_path = "/out/intersection/fdr0.05/5e-08/500000/Kasela_2017_microarray_T-cell_CD8/ieu-a-1162.vcf.bgz"
# eqtl_leads_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.leadpair.tsv"
# eqtl_all_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.all.regions.tsv.gz"
# out_tsv_path = "../../coloc.tsv"

args = commandArgs(trailingOnly=TRUE)
window = as.numeric(args[1])
pcutoff = as.numeric(args[2])
gwas_vcf_path = args[3]
eqtl_leads_path = args[4]
eqtl_all_path = args[5]
out_tsv_path = args[6]

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
out_columns = c("snp", "SNP.PP.H4", "PP.H4.abf", "nsnps", "PP.H3.abf", "PP.H2.abf", "PP.H1.abf", "PP.H0.abf", "coloc_region", "molecular_trait_id", "egene_ensg")
out_df = data.frame(matrix(vector(), 0, length(out_columns), dimnames=list(c(), out_columns)),  stringsAsFactors=F)

################################################################################
# Load lead eqtl pairs
eqtl_leads_df = read.table(eqtl_leads_path, sep='\t', header=T)
eqtl_leads_df = eqtl_leads_df[eqtl_leads_df$chromosome %in% as.character(1:22), ]
eqtl_leads_df$start = eqtl_leads_df$position - (window/2)
eqtl_leads_df$end = eqtl_leads_df$position + (window/2)
eqtl_leads_df[eqtl_leads_df$start<0, "start"] <-0  # negative values with 0
region_lst = as.character(paste0(eqtl_leads_df$chrom, ':', eqtl_leads_df$start, '-', eqtl_leads_df$end))

################################################################################
# Loop over regions
region_i = region_lst[1]
region_i = "1:120860278-121860278"
for (region_i in c("1:120860278-121860278")) {
  print(region_i)

  ################################################################################
  # Load and filter gwas summary statistics
  gwas_vcf <- gwasvcf::query_chrompos_file(chrompos=region_i, vcffile=gwas_vcf_path)
  if(nrow(gwasvcf::vcf_to_tibble(gwas_vcf)) == 0) {next}  # continue if empty gwas
  gwas_vcf <- gwasvcf::query_pval_vcf(vcf=gwas_vcf, pval=pcutoff)  # select associations
  if(nrow(gwasvcf::vcf_to_tibble(gwas_vcf)) == 0) {next}  # continue if empty gwas
  gwas_tbl <- gwas_vcf %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  
  gwas_tbl = gwas_tbl[gwas_tbl$AF<1, ]  # keep MAF<1
  gwas_tbl = gwas_tbl[gwas_tbl$AF>0, ]  # keep MAF>0
  gwas_tbl = gwas_tbl[!duplicated(gwas_tbl$ID), ]  # keep unique RSIDs
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$ES), ]  # keep non-null effect size/Z
  gwas_tbl = gwas_tbl[!is.na(gwas_tbl$SE), ]  # remove non-null se
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
  eqtl_tbl = eqtl_tbl[!duplicated(eqtl_tbl$rsid), ]  # keep non-duplicated RSIDs
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$beta), ]  # keep non-null effect size, z, beta
  eqtl_tbl = eqtl_tbl[!is.na(eqtl_tbl$se), ]  # remove non-null se
  if(nrow(eqtl_tbl) == 0) {next}  # continue if empty gwas
  
  ################################################################################
  # Loop over egenes
  molecular_trait_id_lst = unique(eqtl_tbl$molecular_trait_id)
  molecular_trait_id_i = molecular_trait_id_lst[1]
  
  for (molecular_trait_id_i in molecular_trait_id_lst) {
    
    print(molecular_trait_id_i)
    
    egene_ensg = as.character(eqtl_tbl[head(which(eqtl_tbl$molecular_trait_id==molecular_trait_id_i), 1), "molecular_trait_object_id"])
    # print(egene_ensg)
  
    eqtl_egene_tbl = eqtl_tbl[eqtl_tbl$molecular_trait_id==molecular_trait_id_i, ]
    if(nrow(eqtl_egene_tbl) == 0) {next}  # continue if empty gwas
    
    ################################################################################
    # Coloc
    
    # Keep common SNPs
    rsid_intersection_lst = Reduce(intersect, list(eqtl_egene_tbl$rsid, gwas_tbl$ID))
    if(length(rsid_intersection_lst) == 0) {next}  # continue if empty gwas
    gwas_tbl = gwas_tbl[gwas_tbl$ID %in% rsid_intersection_lst, ]  # keep common rsids
    eqtl_egene_tbl = eqtl_egene_tbl[eqtl_egene_tbl$rsid %in% rsid_intersection_lst, ]  # keep common rsids
    
    # Format for coloc
    type1='quant'
    eqtl_coloc_lst = list(N = (eqtl_egene_tbl$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                          MAF = eqtl_egene_tbl$maf, 
                          beta = eqtl_egene_tbl$beta,
                          varbeta = eqtl_egene_tbl$se^2, 
                          type = type1, 
                          snp = eqtl_egene_tbl$rsid)
    
    
    gwas_coloc_lst = gwas_tbl %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, 
                                        beta = .$ES, varbeta = .$SE^2, type = type1, 
                                        snp = .$ID, z = .$ES / .$SE, 
                                        chr = .$seqnames, pos = .$start, 
                                        id = VariantAnnotation::samples(VariantAnnotation::header(gwas_vcf))[1])}
    
    # print(eqtl_coloc_lst)
    # print(gwas_coloc_lst)
    print("error????????????????")
    coloc_res <- coloc::coloc.abf(eqtl_coloc_lst, gwas_coloc_lst)
    coloc_df = coloc_res$results[, c('snp', 'SNP.PP.H4')]
    coloc_df = coloc_df[coloc_df$SNP.PP.H4>0.2, ]
    
    coloc_df = coloc_df %>% dplyr::mutate(PP.H4.abf=coloc_res$summary[['PP.H4.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(nsnps=coloc_res$summary[['nsnps']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H3.abf=coloc_res$summary[['PP.H3.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H2.abf=coloc_res$summary[['PP.H2.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H1.abf=coloc_res$summary[['PP.H1.abf']])
    coloc_df = coloc_df %>% dplyr::mutate(PP.H0.abf=coloc_res$summary[['PP.H0.abf']])
    
    coloc_df = coloc_df %>% dplyr::mutate(coloc_region = region_i)
    coloc_df = coloc_df %>% dplyr::mutate(molecular_trait_id = molecular_trait_id_i)
    coloc_df = coloc_df %>% dplyr::mutate(egene_ensg = egene_ensg)
    
    coloc_df = dplyr::rename(coloc_df, rsid=snp)
    
    out_df = rbind(out_df, coloc_df)
    # print(coloc_df)
    # print(colnames(coloc_df))
  }
}

# print(out_df)
write.table(out_df, out_tsv_path, quote=F, sep="\t", row.names = F, append=F, col.names = T)
