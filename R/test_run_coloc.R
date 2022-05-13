source('R/coloc_gwasglue_leadpair_all.R')


coloc_window = 500000
eqtl_identifier = "Schmiedel_2018_ge_monocyte_CD16_naive"
eqtlleads_hg38_tsv_path = "results/hg38/eqtlleadpairs/9:4000000-6000000/Schmiedel_2018_ge_monocyte_CD16_naive.leadpairs.tsv"
eqtl_hg38_tsv_gz_path = "ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Schmiedel_2018/ge/Schmiedel_2018_ge_monocyte_CD16_naive.all.tsv.gz"
gwas_identifier = "ieu-a-1008"
gwas_trait_name = "Platelet count"
gwastop_hg38_tsv_path = "results/gwas/5e-08/ieu-a-1008_hg38.tsv"
gwas_hg38_vcf_gz_path = "resources/hg38/gwas.mrcieu.ac.uk/files/ieu-a-1008/ieu-a-1008_hg38_sorted.vcf.gz"
out_tsv_path = "results/coloc/9:4000000-6000000/5e-08/500000/Schmiedel_2018_ge_monocyte_CD16_naive/ieu-a-1008.tsv"

coloc_gwasglue_leadpair_all(
  coloc_window,
  eqtl_identifier,
  eqtlleads_hg38_tsv_path,
  eqtl_hg38_tsv_gz_path,
  gwas_identifier,
  gwas_trait_name,
  gwastop_hg38_tsv_path,
  gwas_hg38_vcf_gz_path,
  out_tsv_path
)
