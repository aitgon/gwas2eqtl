split_into_chunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = floor(n_total/(n_chunks))
  source("R/split_into_batches.R")
  batches = split_into_batches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}
