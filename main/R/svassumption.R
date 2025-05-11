# Suppress package startup messages and load libraries
suppressPackageStartupMessages({
  invisible(library(data.table))
  invisible(library(dplyr))
  invisible(library(coloc))
  invisible(library(hash))
  invisible(library(optparse))
  invisible(library(R.utils))
  invisible(library(stringr))
  invisible(library(parallel))
})



# String concatenation operator
"%&%" <- function(a, b) paste0(a, b)

# Command line arguments.
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", 
              help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--population", type = "character", 
              help = "The subpopulation to create LD matrix (1000 Genomes). See README.md."),
  make_option("--output", type = "character", 
              help = "Output directory."),
  make_option("--data", type = "character", 
              help = "Directory where preprocessed GWAS and QTL data is stored."),
  make_option("--ld_input", type = "character", 
              help = "Directory for LD data."),
  make_option("--threads", type = "character", default = "2", 
              help = "Number of threads to use.")
)

# Parse command line arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list))

# Get list of gene directories in the TWAS folder.
genes_dirs <- list.dirs(opt$data, full.names = TRUE)[-1]

# Run single-variant coloc for each gene.
for (gene_dir in gene_dirs) {
    gene_name <- basename(gene_dir)
    if (length(list.files(member)) == 0) {
        message("SVA: Skipped unprocessed entry (likely missing in QTL lists): ", gene_name)
        next
    }
    saveRDS(sva_QTL(gene_dir), file = gene_dir %&% "/" %&% gene_name %&% "_sva")
}

