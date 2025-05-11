# ----------------------------
# Concatenate and Index VCFs
# ----------------------------

# Concatenate all corrected chromosome files
all_chr_vcf <- paste0(LD_output, "/all_chr.vcf.gz")
if (!file.exists(all_chr_vcf)) {
  input_files <- LD_folder %&% "/ALL.chr*.vcf.gz"
  cmd <- paste(
    "bcftools concat -Oz --threads", opt$threads,
    "-o", shQuote(all_chr_vcf),
    input_files
  )
  message("Running: ", cmd)
  system(cmd)
}

# Create index file for concatenated file.
if (!file.exists(all_chr_vcf %&% ".csi") & !file.exists(all_chr_vcf %&% ".tbi")) {
  message("Indexing VCF...")
  cmd <- paste(
    "tabix -p vcf", shQuote(all_chr_vcf)
  )
  message("Running: ", cmd)
  system(cmd)
}

# Concatenating annotation file.
all_chr_annotation_dir <- paste0(LD_folder, "/annotation")
all_chr_annotation <- paste0(LD_output, "/all_chr_annotation.vcf.gz")
if(!file.exists(all_chr_annotation)) {
  input_files <- all_chr_annotation_dir %&% "/*.vcf.gz"
  cmd <- paste(
    "bcftools concat -Oz --threads", opt$threads,
    "-o", shQuote(all_chr_annotation),
    input_files
  )
  message("Running: ", cmd)
  system(cmd)
}

# Create index file for concatenated file.
if (!file.exists(all_chr_annotation %&% ".csi") & !file.exists(all_chr_annotation %&% ".tbi")) {
  cmd <- paste(
    "tabix -p vcf", shQuote(all_chr_annotation)
  )
  message("Running: ", cmd)
  system(cmd)
}

# Annotate the concatenated VCF
annotated_vcf <- paste0(LD_output, "/all_chr_modified.vcf.gz")
vcf_to_annotate <- LD_output %&% "/all_chr.vcf.gz"
ref_vcf <- all_chr_annotation

# # Annotate with dbSNP if the annotated file doesn't exist
# if (!file.exists(annotated_vcf)) {
#   message("Annotating VCF using renamed files...")
#   cmd <- paste(
#     "bcftools annotate",
#     "-a", shQuote(ref_vcf),
#     "-c ID",
#     "-Oz",                              # Output in bgzipped format
#     "-o", shQuote(annotated_vcf),       # Output file
#     shQuote(vcf_to_annotate),
#     "--threads", opt$threads
#   )
#   message("Running: ", cmd)
#   system(cmd)
# }

# # Index annotated file
# if (!file.exists(annotated_vcf %&% ".csi") & !file.exists(annotated_vcf %&% ".tbi")) {
#   cmd <- paste(
#     "tabix -p vcf", shQuote(annotated_vcf)
#   )
#   message("Running: ", cmd)
#   system(cmd)
# }

# ----------------------------
# Convert VCF to PLINK Format
# ----------------------------