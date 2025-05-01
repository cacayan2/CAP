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

# ----------------------------
# Function: Single-variant coloc for QTL vs GWAS
# ----------------------------
sva_QTL <- function(TWAS_folder, member) {
  member_suffix <- str_replace(paste0("/", member), paste0(TWAS_folder, "/"), "")
  gwascoloc <- readRDS(paste0(member, member_suffix, "_gwascoloc"))
  qtlcoloc  <- readRDS(paste0(member, member_suffix, "_qtlcoloc"))

  sv.res <- coloc.abf(dataset1 = gwascoloc, dataset2 = qtlcoloc)
  return(sv.res)
}

# String concatenation operator
"%&%" <- function(a, b) paste0(a, b)

# ----------------------------
# Command-line argument parsing
# ----------------------------
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", 
              help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--superpop", type = "character", 
              help = "The subpopulation to create LD matrix (1000 Genomes). See README.md."),
  make_option("--output", type = "character", 
              help = "Output directory."),
  make_option("--data", type = "character", 
              help = "Directory where preprocessed GWAS and QTL data is stored."),
  make_option("--lddir", type = "character", 
              help = "Directory for LD data."),
  make_option("--lddownload", type = "character", 
              help = "True/False flag for whether LD data was downloaded using the pipeline."),
  make_option("--threads", type = "character", default = "2", 
              help = "Number of threads to use.")
)

opt <- optparse::parse_args(OptionParser(option_list = option_list))

# Set up paths
LD_folder <- opt$lddir
LD_output <- file.path("LD_output")

# ----------------------------
# Prepare TWAS Data & Run SVA
# ----------------------------

# Remove trailing slash from --data path
TWAS_folder <- str_replace(opt$data, "/", "")

# Get list of gene directories in the TWAS folder
data_members <- list.dirs(TWAS_folder, full.names = TRUE)[-1]

# Run single-variant coloc for each gene
for (member in data_members) {
  member_name <- str_replace(member, paste0(TWAS_folder, "/"), "")
  
  if (length(list.files(member)) == 0) {
    message("SVA: Skipped unprocessed entry (likely missing in QTL lists): ", member_name)
    next
  }

  current_sva <- sva_QTL(TWAS_folder, member)
  saveRDS(current_sva, file = file.path(member, member_name %&% "_sva"))
}

# ----------------------------
# LD Matrix Preprocessing
# ----------------------------

# Create LD output directories if they don’t exist
if (!file.exists(LD_output)) dir.create(LD_output)

# ----------------------------
# Concatenate and Index VCFs
# ----------------------------

# Concatenate all corrected chromosome files
message("Concatenating VCFs...")
all_chr_vcf <- paste0(LD_output, "/all_chr.vcf.gz")
if (!file.exists(all_chr_vcf)) {
  input_files <- LD_folder %&% "/ALL.chr*.vcf.gz"
  cmd <- paste(
    "bcftools concat -Oz --threads", opt$threads,
    "-o", shQuote(all_chr_vcf),
    input_files
  )
  system(cmd)
}

# Create index file for concatenated file.
if (!file.exists(all_chr_vcf %&% ".csi") & !file.exists(all_chr_vcf %&% ".tbi")) {
  message("Indexing VCF...")
  cmd <- paste(
    "tabix -p vcf", shQuote(all_chr_vcf)
  )
  system(cmd)
}

# Mapping for chromosomes 1-22
accessions <- c(
  "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10",
  "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11",
  "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10",
  "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
  "NC_000021.9", "NC_000022.11"
)

chroms <- as.character(1:22)

# Setting directories for mappings
contig_map <- LD_output %&% "/contig_map"
contig_headers <- LD_output %&% "/contig_headers"

# Write to contig_map.txt in required format
lines <- paste0("##contig=<ID=", accessions, ",localAlias=", chroms, ">")
writeLines(lines, contig_headers)

header_lines <- readLines(contig_headers)
contig_lines <- grep("^##contig", header_lines, value = TRUE)

accessions <- str_match(contig_lines, "ID=([^,]+)")[,2]
aliases    <- str_match(contig_lines, "localAlias=([^>]+)")[,2]

contig_map_df <- data.frame(accessions, aliases)
write.table(
  contig_map_df, file = LD_output %&% "/contig_map",
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Set file paths
dbsnp_input <- LD_folder %&% "GCF_000001405.25.gz"
dbsnp_output <- LD_output %&% "/GCF_000001405.25_renamed.gz"

# Construct and run bcftools command to rename chromosomes
cmd <- paste(
  "bcftools annotate",
  "--rename-chrs", shQuote(contig_map),
  "-Oz -o", shQuote(dbsnp_output),
  shQuote(dbsnp_input),
  "--threads", opt$threads
)

if(!file.exists(dbsnp_output)) {
  message("Renaming chromosomes...")
  message("Running: ", cmd)
  system(cmd)
}

# Index renamed file
message("Indexing renamed file...")
if (!file.exists(dbsnp_output %&% ".csi") & !file.exists(dbsnp_output %&% ".tbi")) {
  cmd <- paste(
    "tabix -p vcf", shQuote(dbsnp_output)
  )
  system(cmd)
}

# Annotate the concatenated VCF
annotated_vcf <- paste0(LD_output, "/all_chr_modified.vcf.gz")
dbsnp_vcf <- LD_output %&% "/GCF_000001405.25_renamed.gz" 
vcf_to_annotate <- LD_output %&% "/all_chr.vcf.gz"

# Annotate with dbSNP if the annotated file doesn't exist
if (!file.exists(annotated_vcf)) {
  message("Annotating VCF using renamed files...")
  cmd <- paste(
    "bcftools annotate",
    "-a", shQuote(dbsnp_vcf),              # Columns to annotate (adjust as needed)
    "-c CHROM,POS,ID",
    "-Oz",                              # Output in bgzipped format
    "-o", shQuote(annotated_vcf),       # Output file
    shQuote(vcf_to_annotate),
    "--threads", opt$threads
  )
  system(cmd)
}

# Index annotated file
message("Indexing annotated file...")
if (!file.exists(annotated_vcf %&% ".csi") & !file.exists(annotated_vcf %&% ".tbi")) {
  cmd <- paste(
    "tabix -p vcf", shQuote(annotated_vcf)
  )
  system(cmd)
}

# ----------------------------
# Convert VCF to PLINK Format
# ----------------------------

message("Converting VCF to PLINK...")
plink_prefix <- file.path(LD_output, "all_chr_modified")
if (!file.exists(plink_prefix %&% ".bed") ||
    !file.exists(plink_prefix %&% ".bim") ||
    !file.exists(plink_prefix %&% ".fam")) {
  plink_cmd <- paste(
    "plink2 --vcf", shQuote(all_chr_vcf),
    "--make-bed --out", shQuote(plink_prefix),
    "--allow-extra-chr --max-alleles 2", "--threads", opt$threads
  )
  system(plink_cmd)
}

# Load sample information
message("Loading sample information...")
sample_info <- fread(LD_folder %&% "integrated_call_samples_v3.20200731.ALL.ped")

# Reformat sample data to include PLINK .psam fields and superpopulation information
message("Reformatting sample information...")
sample_info <- sample_info %>% 
  mutate(
    FID = 0,
    IID = `Individual ID`,
    PAT = 0,
    MAT = 0,
    SEX = case_when(
      Gender == "male" ~ 1,
      Gender == "female" ~ 2,
      TRUE ~ 0
    ),
    Population = Population,
    SuperPop = case_when(
      Population %in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB") ~ "AFR",
      Population %in% c("MXL", "PUR", "CLM", "PEL") ~ "AMR",
      Population %in% c("CHB", "JPT", "CHS", "CDX", "KHV") ~ "EAS",
      Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS") ~ "EUR",
      Population %in% c("GIH", "PJL", "BEB", "STU", "ITU") ~ "SAS"
    ),
    PHENOTYPE = -9
  ) %>%
  select(FID, IID, PAT, MAT, SEX, SuperPop, Population, PHENOTYPE)

# Write formatted sample info to .psam file
fwrite(sample_info, LD_output %&% "/all_chr.psam", col.names = TRUE, sep = "\t")

# Load .psam file to get list of individuals
pops <- fread(LD_output %&% "/all_chr.psam")

# Subset individuals based on target superpopulation or population
if (opt$superpop %in% unique(sample_info$Population)) {
  superpop <- filter(pops, Population == opt$superpop)
} else if (opt$superpop == "ALL") {
  superpop <- pops
} else {
  superpop <- filter(pops, SuperPop == opt$superpop)
}

# Write FID/IID list for selected individuals
superpoplist <- mutate(superpop, FID = 0) %>% select(IID)
fwrite(superpoplist, paste0(LD_output, "/superpop_list"), col.names = FALSE, sep = "\t")

# Generate LD matrices and coloc inputs if running for QTL-based analyses
message("Generating LD matrices and coloc inputs...")
if (opt$process == "pqtl" | opt$process == "eqtl") {
  for (member in data_members) {
    # Skip if folder is empty
    if (length(list.files(member)) == 0) next

    # Load coloc input files
    gwascoloc <- readRDS(paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascoloc"))
    qtlcoloc  <- readRDS(paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcoloc"))
    
    # Generate SNP list using annotated .bim file.
    split_pos <- do.call(rbind, strsplit(qtlcoloc$snp, ":"))
    write.table(split_pos, file = paste0(member, "/positions_file"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    stopifnot(file.exists(paste0(member, "/positions_file")))

    # Annotate using bcftools to get rsIDs
    positions_file <- paste0(member, "/positions_file")
    dbsnp_vcf <- paste0(LD_output, "/GCF_000001405.25_renamed.gz")
    all_chr_file <- paste0(LD_output, "/all_chr_modified.vcf.gz")
    output_vcf <- tempfile(fileext = ".vcf.gz")

    cmd <- paste(
      "bcftools query",
      "-R", shQuote(positions_file),
      "-f '%CHROM\t%POS\t%ID\n'",
      shQuote(dbsnp_vcf),
      " > ", member %&% "/matches.txt"
    )
    
    message("Running: ", cmd)
    system(cmd)
    

    positions <- read.table(positions_file, header = FALSE, col.names = c("CHROM", "POS"))

    vcf_matches <- read.table(member %&% "/matches.txt", header = FALSE, col.names = c("CHROM", "POS", "ID"))

    final_output <- merge(vcf_matches, positions, by = c("CHROM", "POS"), all.x = TRUE)

    write.table(final_output$ID, file = paste0(member, "/coloc-snplist"), col.names = FALSE, row.names = FALSE, quote = FALSE)

    # Generate genotype matrix using PLINK
    system("plink2 --bfile " %&% LD_output %&% "/all_chr_modified --extract " %&% member %&% 
           "/coloc-snplist --keep " %&% LD_output %&% 
           "/superpop_list --recode A --make-just-bim --out " %&% 
           member %&% "/superpop_coloc-region --allow-extra-chr --maf 0.01 --threads " %&% opt$threads)

    # Read PLINK genotype data
    geno <- fread(member %&% "/superpop_coloc-region.raw")
    genomat <- as.matrix(geno[, 7:ncol(geno)])
    snplist <- substr(colnames(genomat), 1, nchar(colnames(genomat)) - 2)
    colnames(genomat) <- snplist

    # Read BIM file
    bim <- fread(member %&% "/superpop_coloc-region.bim")

    # Allele complement map
    bases <- hash()
    bases[["A"]] <- "T"
    bases[["T"]] <- "A"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"

    # Match SNPs and handle allele matching / flipping / complementing
    matchsnps <- readRDS(paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_matchsnps"))
    susiesnpsbim <- inner_join(matchsnps, bim, by = c("ID" = "V2"))

    # Direct match
    matchbim <- filter(susiesnpsbim, A1.x == V5 & REF == V6)

    # Complement match
    b <- filter(susiesnpsbim, !(REF == "A" & ALT == "T") & !(REF == "T" & ALT == "A") &
                                  !(REF == "C" & ALT == "G") & !(REF == "G" & ALT == "C"))
    compmatchbim <- filter(b, A1.x == values(bases, keys = V5) & REF == values(bases, keys = V6))

    # Flipped alleles
    flippedbim <- filter(b, A1.x == V6 & REF == V5)
    compflippedbim <- filter(b, A1.x == values(bases, keys = V6) & REF == values(bases, keys = V5))

    # Flip betas if needed
    if (nrow(flippedbim) > 0) flippedbim <- mutate(flippedbim, BETA.x = -BETA.x, BETA.y = -BETA.y)
    if (nrow(compflippedbim) > 0) compflippedbim <- mutate(compflippedbim, BETA.x = -BETA.x, BETA.y = -BETA.y)

    # Combine matched SNPs
    matchbimsnps <- bind_rows(matchbim, compmatchbim, flippedbim, compflippedbim)

    # Create LD matrix
    x <- genomat[, colnames(genomat) %in% matchbimsnps$ID]
    R <- cor(x)

    # Create coloc objects with LD for SuSiE
    susiesnps <- filter(matchbimsnps, ID %in% colnames(x))
    gwascolocsusie <- list(
      beta = susiesnps$BETA.y,
      varbeta = (susiesnps$SE.y)^2,
      snp = susiesnps$ID,
      position = susiesnps$POS,
      type = "cc",
      LD = R,
      N = 100000
    )

    qtlcolocsusie <- list(
      beta = setNames(susiesnps$BETA.x, susiesnps$ID),
      varbeta = (susiesnps$SE.x)^2,
      snp = susiesnps$ID,
      position = susiesnps$POS,
      type = "quant",
      N = susiesnps$OBS_CT[1],
      MAF = susiesnps$A1_FREQ,
      sdY = 1,
      LD = R
    )

    # Save coloc and LD results
    saveRDS(gwascolocsusie, file = paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascolocsusie"))
    saveRDS(qtlcolocsusie,  file = paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcolocsusie"))
    saveRDS(R,              file = paste0(member, str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_R"))
  }
}
