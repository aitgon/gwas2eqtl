source("R/coloc_gwasglue_leadpair_one.R")

eqtl_identifier = "Kasela_2017_microarray_T-cell_CD8"
gwas_identifier = "ieu-a-1162"
lead_egene = "ILMN_2204664"
lead_chrom = 1
lead_pos_hg38 = 121108581
coloc_window = 500000
eqtl_tsv_gz_path = "ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.all.tsv.gz"
gwas_hg38_vcf_gz_path = "resources/hg38/gwas.mrcieu.ac.uk/files/ieu-a-1162/ieu-a-1162_hg38_sorted.vcf.gz"

coloc_df = coloc_gwasglue_leadpair_one(coloc_window, eqtl_identifier, lead_egene, lead_chrom, lead_pos_hg38, eqtl_tsv_gz_path, gwas_identifier, gwas_hg38_vcf_gz_path)
