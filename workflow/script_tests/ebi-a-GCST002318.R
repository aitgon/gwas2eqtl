library(dplyr)

gwas_identifier = "ebi-a-GCST002318"

df = ieugwasr::tophits(gwas_identifier, clump = 1, r2=0.1, kb=1000, pop="EUR")
df = df[order(df$chr, df$position), ]

rsid = "rs28411352"
df0 = df[df$rsid==rsid,]


granges_df = df
granges_df$start = df$position
granges_df$end = df$position
GenomicRanges::makeGRangesFromDataFrame(granges_df, keep.extra.columns = T)

df <-df[order(df$chr, df$position),]
df

eqtl_tbl = seqminer::tabix.read.table(tabixFile = eqtl_all_path, tabixRange = region_i, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

gwasvcf::query_pval_file(vcf="/home/gonzalez/Software/process/hg38/gwas.mrcieu.ac.uk/files/ebi-a-GCST002318/ebi-a-GCST002318.vcf.bgz", pval=5e-8)

gwasvcf::query_chrompos_file(chrompos=region_i, vcffile=gwas_vcf_path)