

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
  sv.sens <- sensitivity(sv.res, "H4 > 0.5")

  return(c(sv.res, sv.sens))
}

# Define command-line options.
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--superpop", type = "character", help = "The subpopulation to create LD matrix (generated from 1000 genomes project). See README.md for more options.")
)

# Obtain command line arguments.
opt <- parse_args(OptionParser(option_list = option_list))
if(length(opt) != 3) {
  print("Please pass an appropriate number of arguments to the command line.", quote = FALSE)
  stop()
}

# Set folder for TWAS data. 
TWAS_folder = ""

if(opt$process == "eQTL") {
  TWAS_folder = paste0(getwd(), "/processed/eQTL_results")
} else {
  TWAS_folder = paste0(getwd(), "/processed/pQTL_results")
}

message("Running coloc assuming single variant assumption...")
# Running coloc using single variant assumption
data_members <- list.dirs(TWAS_folder)
data_members <- data_members[-1]

sv.analysis <- new.env()
for(member in data_members) {
  if(length(list.files(member)) == 0) {
    message(paste0("SVA: The following result not processed - likely not an entry in QTL lists: ", str_replace(member, paste0(TWAS_folder, "/"), "")))
    next
  }
  sv.analysis[[str_replace(member, paste0(TWAS_folder, "/"), "")]] <- sva_QTL(TWAS_folder, member)
}

print("Setting up linkage disequilibrium data...")
