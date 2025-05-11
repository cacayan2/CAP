sva_QTL <- function(gene_dir) {
  """
  Function used to perform single variant analysis on data. 
  """
  gene_name <- basename(gene_dir)
  gwascoloc <- readRDS(gene_dir %&% "/" %&% gene_name %&% "_gwascoloc")
  qtlcoloc  <- readRDS(gene_dir %&% "/" %&% gene_name %&% "_qtlcoloc")

  sv.res <- coloc.abf(dataset1 = gwascoloc, dataset2 = qtlcoloc)
  return(sv.res)
}