source("R/coloc_gwasglue_leadpair_one.R")

eqtl_identifier = "Schmiedel_2018_ge_Treg_naive"
gwas_identifier = "ebi-a-GCST004603"
lead_egene = "ENSG00000106686"
lead_chrom = 9
lead_pos_hg38 = 4676745
coloc_window = 500000
eqtl_tsv_gz_path = "ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Schmiedel_2018/ge/Schmiedel_2018_ge_Treg_naive.all.tsv.gz"
gwas_hg38_vcf_gz_path = "resources/hg38/gwas.mrcieu.ac.uk/files/ebi-a-GCST004603/ebi-a-GCST004603_hg38_sorted.vcf.gz"

coloc_df = coloc_gwasglue_leadpair_one(coloc_window, eqtl_identifier, lead_egene, lead_chrom, lead_pos_hg38, eqtl_tsv_gz_path, gwas_identifier, gwas_hg38_vcf_gz_path)
