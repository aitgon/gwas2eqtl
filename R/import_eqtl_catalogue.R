library(dplyr)
# In eQTL Catalogue, **variants with multiple rsids are split over multiple rows** in the summary statistics files.
# Thus, we first want to retain only one unique record per variant. 
# To simplify colocalisation analysis, we also want to exclude multi-allelic variants. 
# The following function imports summary statistics from a tabix-index TSV file and performs necessary filtering.
import_eqtl_catalogue <- function(eqtl_sumstats_path, region, molecular_trait_id, verbose = TRUE){
  # print(cat(eqtl_sumstats_path, region, molecular_trait_id, column_names))
  #Fetch summary statistics with seqminer
  col_names = colnames(read.table(eqtl_sumstats_path, nrows=1, header = T, sep = "\t"))
  eqtl_sumstats_hg38_tbl = seqminer::tabix.read.table(tabixFile = eqtl_sumstats_path, tabixRange = region, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  colnames(eqtl_sumstats_hg38_tbl) = col_names
  cat("eqtl_sumstats_hg38_tbl 1: ", dim(eqtl_sumstats_hg38_tbl), "\n")

  # select lead egene
  eqtl_sumstats_hg38_tbl = eqtl_sumstats_hg38_tbl[eqtl_sumstats_hg38_tbl$molecular_trait_id==molecular_trait_id, ]
  #cat("eqtl_sumstats_hg38_tbl 2: ", dim(eqtl_sumstats_hg38_tbl), "\n")

  eqtl_sumstats_hg38_tbl$rsid = stringr::str_replace_all(eqtl_sumstats_hg38_tbl$rsid, "[\r\n]" , "")
  #cat("eqtl_sumstats_hg38_tbl 3: ", dim(eqtl_sumstats_hg38_tbl), "\n")
  
  if (is.null(eqtl_sumstats_hg38_tbl)) {
    return(NULL)
  }
  eqtl_sumstats_hg38_tbl = eqtl_sumstats_hg38_tbl %>% dplyr::filter(molecular_trait_id == molecular_trait_id)

  # remove variants with unknown SE
  eqtl_sumstats_hg38_tbl = eqtl_sumstats_hg38_tbl %>% dplyr::filter(!is.nan(se))

  ##Remove rsid duplicates and multi-allelic variant
  #  eqtl_sumstats_hg38_tbl = dplyr::select(eqtl_sumstats_hg38_tbl, everything()) %>%
  #  dplyr::distinct() %>% #rsid duplicates
  #  dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
  #  dplyr::group_by(id) %>%
  #  dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>%
  #  dplyr::filter(row_count == 1) %>% #Multialllics
  #  dplyr::filter(!is.nan(se)) %>% # remove variants with unknown SE
  #  dplyr::filter(!is.na(se)) %>%
  #  dplyr::select(-row_count) # remove added row_count columns
  return(eqtl_sumstats_hg38_tbl)
}
