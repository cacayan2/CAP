

# Loading appropriate libraries. 
suppressPackageStartupMessages(invisible(library(data.table)))
suppressPackageStartupMessages(invisible(library(dplyr)))
suppressPackageStartupMessages(invisible(library(coloc)))
suppressPackageStartupMessages(invisible(library(hash)))
suppressPackageStartupMessages(invisible(library(optparse)))
suppressPackageStartupMessages(invisible(library(R.utils)))
suppressPackageStartupMessages(invisible(library(stringr)))
suppressPackageStartupMessages(invisible(library(hash)))

sva_QTL <- function(TWAS_folder, member) {
  gwascoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascoloc")))
  qtlcoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcoloc")))

  sv.res <- coloc.abf(dataset1 = gwascoloc, dataset2 = qtlcoloc)
}

# Define command-line options.
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--superpop", type = "character", help = "The subpopulation to create LD matrix (generated from 1000 genomes project). See README.md for more options."),
  make_option("--output", type = "character", help = "Output directory.")
)

# Obtain command line arguments.
opt <- parse_args(OptionParser(option_list = option_list))
if(length(opt) != 4) {
  print("Please pass an appropriate number of arguments to the command line.", quote = FALSE)
  stop()
}

# Set folder for TWAS data. 
TWAS_folder = paste0(getwd(), paste0("/", opt$output))


message("Running coloc assuming single variant assumption...")
# Running coloc using single variant assumption
data_members <- list.dirs(TWAS_folder)
data_members <- data_members[-1]
for(member in data_members) {
  member_name = str_replace(member, paste0(TWAS_folder, "/"), "")
  if(length(list.files(member)) == 0) {
    message(paste0("SVA: The following result not processed - likely not an entry in QTL lists: ", member_name))
    next
  }
  current_sva <- sva_QTL(TWAS_folder, member)
  saveRDS(current_sva, file = file.path(member, paste0(member_name, "_sva"))) 
}

message("Setting up linkage disequilibrium data...")
