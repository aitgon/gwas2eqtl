coloc_gwasglue_leadpair_all <- function(coloc_window, eqtl_identifier, eqtlleads_hg38_tsv_path, eqtl_hg38_tsv_gz_path, gwas_identifier, gwas_trait_name, gwastop_hg38_tsv_path, gwas_hg38_vcf_gz_path, out_tsv_path) {
  
  ################################################################################
  #
  # libraries
  #
  ################################################################################
  
  suppressPackageStartupMessages(suppressWarnings({
    library(gwasglue)
    library(dplyr)
    library(gassocplot)
    library(coloc)
    library(readr)
    library(GenomicRanges)
  }))
  source("R/import_eqtl_catalogue.R")
  source("R/split_into_chunks.R")
  source("R/coloc_gwasglue_leadpair_one.R")
  
  ################################################################################
  #
  # script
  #
  ################################################################################
  
  #Delete file if it exists
  if (file.exists(out_tsv_path)) {file.remove(out_tsv_path)}
  
  # Create output dir
  if (!dir.exists(dirname(out_tsv_path))) dir.create(dirname(out_tsv_path))

  # Read molecular_trait_id variant pairs and convert it to granges
  eqtlleads_tbl <- readr::read_tsv(eqtlleads_hg38_tsv_path)
  colnames(eqtlleads_tbl) = c("egene", "variant", "chr", "position")
  
  # eqtleands has no rows
  if (dim(eqtlleads_tbl)[[1]]==0) {
    cat("No significant lead eQTLs\n")
    cat(NULL, file=out_tsv_path)
    return(NULL)
  }
  
  eqtlleads_tbl$start = eqtlleads_tbl$position - coloc_window
  eqtlleads_tbl$start[eqtlleads_tbl$start<0] <- 0
  eqtlleads_tbl$end = eqtlleads_tbl$position + coloc_window
  eqtlleads_grngs = GenomicRanges::makeGRangesFromDataFrame(eqtlleads_tbl, keep.extra.columns = T)

  # top gwasglue and convert them to granges
  # gwastop_tbl <- ieugwasr::tophits(gwas_identifier, pval = gwas_pval)
  if (file.info(gwastop_hg38_tsv_path)$size==0) {
    cat("No significant GWAS signals at p-value\n")
    cat(NULL, file=out_tsv_path)
    return(NULL)
  }

  gwastop_tbl <- readr::read_tsv(gwastop_hg38_tsv_path)
  gwastop_tbl = gwastop_tbl[!is.na(gwastop_tbl$n),]  # remove n NAs
  if (nrow(gwastop_tbl) == 0) { # if non top gwas, then write empty file
    cat("No significant GWAS signals at p-value\n")
    cat(NULL, file=out_tsv_path)
    return(NULL)
  }

  gwastop_tbl <- gwastop_tbl %>% arrange(p)
  gwastop_tbl$start = gwastop_tbl$position - coloc_window
  gwastop_tbl$start[gwastop_tbl$start<0] <- 0
  gwastop_tbl$end = gwastop_tbl$position + coloc_window
  gwastop_grngs = GenomicRanges::makeGRangesFromDataFrame(gwastop_tbl, keep.extra.columns = T)
  
  # Find overlapping ranges
  overlap_tbl = dplyr::as_tibble(GenomicRanges::findOverlaps(eqtlleads_grngs, gwastop_grngs, select="all", minoverlap=1))
  print(overlap_tbl)
  selected_pairs_tbl = dplyr::as_tibble(eqtlleads_grngs[unique(overlap_tbl$queryHits)])
  print(selected_pairs_tbl)
  selected_gwastop_tbl = dplyr::as_tibble(gwastop_grngs[unique(overlap_tbl$subjectHits)])
  #print(selected_gwastop_tbl[c('seqnames', 'position', 'rsid', 'start', 'end')])

  if (nrow(selected_pairs_tbl) == 0){  # if empty overlapping, write empty file
    cat("No overlapping between eQTL leads and GWAS signals", "\n")
    cat(NULL, file=out_tsv_path)
    return(NULL)
  }
  
  #colnames(selected_pairs_tbl) = c("chromosome", "start", "end", "width", "strand", "egene", "variant", "position")

  cat("nrow selected_pairs_tbl", nrow(selected_pairs_tbl))
  for(i in 1:nrow(selected_pairs_tbl)) {
    #cat("\n\n")
    #cat("i-th eQTL: ", i, ". Total nb of eqtls: ",  nrow(selected_pairs_tbl), "\n")
    lead_egene <- selected_pairs_tbl[[i, "egene"]]
    lead_chrom <- selected_pairs_tbl[[i, "seqnames"]]
    lead_pos_hg38 <- selected_pairs_tbl[[i, "position"]]
    #cat(paste0("Lead egene: ", lead_egene, ", lead chrom: ", lead_chrom, ", lead position hg38: ", lead_pos_hg38, "\n"))

    out_df <- coloc_gwasglue_leadpair_one(coloc_window, eqtl_identifier, lead_egene, lead_chrom, lead_pos_hg38, eqtl_hg38_tsv_gz_path, gwas_identifier, gwas_trait_name, gwas_hg38_vcf_gz_path)
    if (!is.null(out_df)) {
      if (file.exists(out_tsv_path)) { # output file exists
        write.table(out_df, out_tsv_path, quote=F, sep="\t", row.names = F, append=T, col.names = F)
      } else {
        write.table(out_df, out_tsv_path, quote=F, sep="\t", row.names = F, append=T, col.names = T)
      }
    } 
  }
  if (!file.exists(out_tsv_path)) {
    cat("No colocalisation results found", "\n")
    cat(NULL, file=out_tsv_path)
  }
}
