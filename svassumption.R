suppressPackageStartupMessages(invisible(library(data.table)))
suppressPackageStartupMessages(invisible(library(dplyr)))
suppressPackageStartupMessages(invisible(library(coloc)))
suppressPackageStartupMessages(invisible(library(hash)))
suppressPackageStartupMessages(invisible(library(optparse)))
suppressPackageStartupMessages(invisible(library(R.utils)))
suppressPackageStartupMessages(invisible(library(stringr)))
suppressPackageStartupMessages(invisible(library(parallel)))
suppressPackageStartupMessages(invisible(library(parallelly)))

# Function that loads the datasets for single variant assumption. 
sva_QTL <- function(TWAS_folder, member) {
  gwascoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascoloc")))
  qtlcoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcoloc")))

  sv.res <- coloc.abf(dataset1 = gwascoloc, dataset2 = qtlcoloc)
}

# Define string concatenation function.
"%&%" <- function(a,b) paste(a,b,sep="")

# Define command-line options.
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--superpop", type = "character", help = "The subpopulation to create LD matrix (generated from 1000 genomes project). See README.md for more options."),
  make_option("--output", type = "character", help = "Output directory."),
  make_option("--data", type = "character", help = "Directory where data for gwas and eqtl from preprocessing is stored."),
  make_option("--lddir", type = "character", help = "Directory for LD data."),
  make_option("--lddownload", type = "character", help = "True/False whether or not LD data was downloaded using the pipeline."),
  make_option("--threads", type = "character", default = "2", help = "Number of threads to use.")
)

# Obtain command line arguments.
opt <- parse_args(OptionParser(option_list = option_list))
LD_folder = opt$lddir
LD_output = "LD_output/"

# Set folder for TWAS data. 
TWAS_folder = str_replace(opt$data, "/", "")

# Create list of genes in the TWAS folder. 
data_members <- list.dirs(TWAS_folder)
data_members <- data_members[-1]

# Run SVA for each gene. 
for(member in data_members) {
  member_name = str_replace(member, TWAS_folder %&% "/", "")
  if(length(list.files(member)) == 0) {
    message(paste0("SVA: The following result not processed - likely not an entry in QTL lists: ", member_name))
    next
  }
  current_sva <- sva_QTL(TWAS_folder, member)
  saveRDS(current_sva, file = file.path(member, member_name %&% "_sva")) 
}

# Create folder to contain output for linkage disequilibrium processing. 
if(!file.exists(LD_output)) {
  dir.create(LD_output)
}

if(!dir.exists(LD_output %&% "corrected")) {
  dir.create(LD_output %&% "corrected")
}

# In R with %&%
files = list.files("vcf_folder", pattern = "\\.vcf\\.gz$", full.names = TRUE)

# Now, scale bcftools usage based on the threads
mclapply(list.files(LD_folder, full.names = FALSE)[-1], function(f) {
  corrected_file <- LD_output %&% "corrected/" %&% f %&% "_corrected.vcf.gz"
  
  if (!file.exists(corrected_file)) {
    message("Processing:", f, "\n")  # Debugging print
    
    # Scale the number of threads used by bcftools
    # If opt$threads is less than the available cores, it will use opt$threads
    cmd1 <- "bcftools view -m2 -M2 -i \"MAF>0.01\" -Oz -o \"" %&%
            corrected_file %&% "\" " %&%
            LD_folder %&% f %&% " --threads " %&% opt$threads  # Set to opt$threads
    system(cmd1)
  }
}, mc.cores = as.integer(opt$threads))

# Concatenating all corrected .vcf files to a single larger files. 
if (!file.exists(LD_output %&% "all_chr.vcf")) {
  input_files <- Sys.glob(LD_output %&% "corrected/ALL.chr*_corrected.vcf.gz")
  input_str <- paste(input_files, collapse = " ")
  
  cmd <- "bcftools concat -Ov --threads " %&%
         opt$threads %&% " -o " %&% LD_output %&%
         "all_chr.vcf " %&% input_str
  
  system(cmd)
}

# Creating .bed., .bim., .fam. files for all individuals in all_phase3_combined.vcf.gz. 
if(!file.exists(LD_output %&% "all_chr.bed") | !file.exists(LD_output %&% "all_chr.bim") | !file.exists(LD_output %&% "all_chr.fam")) {
  system("plink2 --vcf " %&% LD_output %&% "all_phase3_combined.vcf --make-bed --out " %&% LD_output %&% "all_chr --allow-extra-chr --threads " %&% opt$threads)
}


# Obtaining sample info. 
sample_info <- fread(LD_folder %&% "20130606_sample_info.txt")

# Reformatting sample info data, adding superpopulation information. 
sample_info <- sample_info %>% 
  mutate(
    FID = 0, 
    IID = Sample, 
    PAT = 0, 
    MAT = 0,
    SEX = case_when(
      Gender == "male" ~ 1,
      Gender == "female" ~ 2,
      TRUE ~ 0
    ),
    Population = Population,
    SuperPop = case_when(
      Population == "YRI" ~ "AFR",
      Population == "LWK" ~ "AFR",
      Population == "GWD" ~ "AFR",
      Population == "MSL" ~ "AFR",
      Population == "ESN" ~ "AFR",
      Population == "ASW" ~ "AFR",
      Population == "ACB" ~ "AFR",
      Population == "MXL" ~ "AMR",
      Population == "PUR" ~ "AMR",
      Population == "CLM" ~ "AMR",
      Population == "PEL" ~ "AMR",
      Population == "CHB" ~ "EAS",
      Population == "JPT" ~ "EAS",
      Population == "CHS" ~ "EAS",
      Population == "CDX" ~ "EAS",
      Population == "KHV" ~ "EAS",
      Population == "CEU" ~ "EUR",
      Population == "TSI" ~ "EUR",
      Population == "FIN" ~ "EUR",
      Population == "GBR" ~ "EUR",
      Population == "IBS" ~ "EUR",
      Population == "GIH" ~ "SAS",
      Population == "PJL" ~ "SAS",
      Population == "BEB" ~ "SAS",
      Population == "STU" ~ "SAS",
      Population == "ITU" ~ "SAS"
    ),
    PHENOTYPE = -9
  ) %>%
  select(FID, IID, PAT, MAT, SEX, SuperPop, Population, PHENOTYPE)
fwrite(sample_info, LD_output %&% "all_chr.psam", col.names = TRUE, sep = "\t")

# Obtaining list of all individuals. 
pops <- fread(LD_output %&% "all_chr.psam")

# Filtering given the population of interest. 
if (opt$superpop %in% unique(sample_info$Population)) {
  superpop <- filter(pops, Population == opt$superpop)
} else {
  superpop <- filter(pops, SuperPop == opt$superpop)
}

# Creating list of all individuals in the population of interest.
superpoplist <- mutate(superpop, FID=0) %>% select(FID, IID)
fwrite(superpoplist, paste0(LD_output, "superpop_list"), col.names = FALSE, sep = "\t")

# Generating datasets for multivariant assumption alongside correlation matrix. 
if(opt$process == "pqtl") {
  for (member in data_members) {
    if(length(list.files(member)) == 0) {
    next
    }
    gwascoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascoloc")))
    qtlcoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcoloc")))
    fwrite(data.frame(qtlcoloc$snp), file = paste0(member, "/coloc-snplist"), col.names = FALSE)
    print(member)
    system("plink2 --bfile " %&% LD_output %&% "all_chr --extract \"" %&% member %&% "/coloc-snplist\" --keep " %&% LD_output %&% "superpop_list --maf 0.01 --recode A --make-just-bim --out \""  %&% member %&% "/superpop_coloc-region\" --allow-extra-chr --threads " %&% opt$threads)
    geno <- fread(member %&% "/superpop_coloc-region.raw")
    genomat <- as.matrix(geno[,7:length(geno)])
    snplist <- colnames(genomat)
    snplist <- substr(snplist, 1, nchar(snplist)-2)
    colnames(genomat) <- snplist
    bim <- fread(member %&% "superpop_coloc-region.bim")
    
    joined <- inner_join(gwascoloc, qtlcoloc, by = c("ID" = "RS"))

    bases = hash()
    bases[["A"]] <- "T"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"
    bases[["T"]] <- "A"
    
    match <- joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

    a <- joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
                 !(joined$REF=="G" & joined$ALT=="C") ]

    compmatch <- a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]
    flipped <- a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
    compflipped <- a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

    if(dim(flipped)[1] > 0){
      flipped = mutate(flipped,BETA.y = -1*BETA.y)
    }
    if(dim(compflipped)[1] > 0){
      compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
    }
    matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(chr_pos)
    
    susiesnpsbim <- inner_join(matchsnps, bim, by = c("ID" = "V2"))
    
    matchbim <- susiesnpsbim[(susiesnpsbim$A1.x==susiesnpsbim$V5 & susiesnpsbim$REF == susiesnpsbim$V6),]

    b <- susiesnpsbim[!(susiesnpsbim$REF=="A" & susiesnpsbim$ALT=="T") & !(susiesnpsbim$REF=="T" & susiesnpsbim$ALT=="A") & !(susiesnpsbim$REF=="C" & susiesnpsbim$ALT=="G") &
                        !(susiesnpsbim$REF=="G" & susiesnpsbim$ALT=="C") ]

    compmatchbim <- b[(b$A1.x==values(bases,keys=b$V5) & b$REF == values(bases,keys=b$V6)),]
    flippedbim <- b[(b$A1.x==b$V6 & b$REF == b$V5),]
    compflippedbim <- b[(b$A1.x==values(bases,keys=b$V6) & b$REF == values(bases,keys=b$V5)),]

    if(dim(flippedbim)[1] > 0){
      flippedbim <- mutate(flippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
    }
    if(dim(compflippedbim)[1] > 0){
      compflippedbim <- mutate(compflippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
    }

    matchbimsnps <- rbind(matchbim, compmatchbim, flippedbim, compflippedbim) %>% arrange(POS)

    snplist <- matchbimsnps$ID

    x <- genomat[,colnames(genomat) %in% snplist]
    R <- cor(x)

    susiesnps <- filter(matchbimsnps, ID %in% snplist)
    gwascolocsusie <- list("beta" = susiesnps$BETA.y, "varbeta" = (susiesnps$SE.y)^2, "snp" = susiesnps$ID, "position" = susiesnps$POS, "type" = "cc", "LD" = R, "N" = 100000)
    qtlcolocsusie <- list("beta" = setNames(susiesnps$BETA.x, susiesnps$ID), "varbeta" = (susiesnps$SE.x)^2, "snp" = susiesnps$ID, "position" = susiesnps$POS, "type" = "quant", "N" = susiesnps$OBS_CT[1], "MAF" = susiesnps$A1_FREQ, "sdY"=1, LD = R)
    saveRDS(gwascolocsusie, file = file.path(parent_dir, target_gene, paste0(target_gene, "_gwascolocsusie")))
    saveRDS(qtlcolocsusie, file = file.path(parent_dir, target_gene, paste0(target_gene, "_qtlcolocsusie")))
    saveRDS(R, file = file.path(parent_dir, target_gene, paste0(target_gene, "_LDcorr")))
  } 
} else {
  for (member in data_members) {
    gwascoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_gwascoloc")))
    qtlcoloc <- readRDS(paste0(member, paste0(str_replace(paste0("/", member), paste0(TWAS_folder, "/"), ""), "_qtlcoloc")))
    fwrite(data.frame(qtlcoloc$snp), file = file.path(member, "coloc-snplist"), col.names = FALSE)
    bcf_command = "bcftools query -f '%CHROM\t%pos\t%ID\n' -R " %&% member %&% "coloc-snplist" %&% LD_output %&% "all_phase3_combined.vcf.gz > "%&% member %&%"snp_pos_to_rsid.tsv --threads " %&% opt$threads
    system(bcf_command)
    eqtl <- fread(member %&%"snp_pos_to_rsid.tsv")
    snpmap <- fread(member %&%"snp_pos_to_rsid.tsv", col.names = c("CHR", "POS", "rsid"))
    eqtl_mapped <- inner_join("snp_pos_to_rsid.tsv", snpmap, by=c("CHR", "POS"))
    fwrite(data.frame(eqtl_mapped$rsid), file = file.path(member, "coloc-snplist"), col.names = FALSE)
    system("plink --bfile " %&% LD_output %&% "all_phase3 --extract " %&% member %&% "coloc-snplist --keep" %&% LD_output %&% "superpop_list --maf 0.01 --make-just-bim --out " %&% member %&% "superpop_coloc-region --allow-extra-chr --threads " %&% opt$threads)
    geno <- fread(member %&% "superpop_coloc-region.raw")
    genomat <- as.matrix(geno[,7:length(geno)])
    snplist <- colnames(genomat)
    snplist <- substr(snplist, 1, nchar(snplist)-2)
    colnames(genomat) <- snplist
    bim <- fread(member %&% "superpop_coloc-region.bim")
    
    joined <- inner_join(gwascoloc, qtlcoloc, by = c("ID" = "RS"))

    bases = hash()
    bases[["A"]] <- "T"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"
    bases[["T"]] <- "A"
    
    match <- joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

    a <- joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
                 !(joined$REF=="G" & joined$ALT=="C") ]

    compmatch <- a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]
    flipped <- a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
    compflipped <- a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

    if(dim(flipped)[1] > 0){
      flipped = mutate(flipped,BETA.y = -1*BETA.y)
    }
    if(dim(compflipped)[1] > 0){
      compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
    }
    matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(chr_pos)
    
    susiesnpsbim <- inner_join(matchsnps, bim, by = c("ID" = "V2"))
    
    matchbim <- susiesnpsbim[(susiesnpsbim$A1.x==susiesnpsbim$V5 & susiesnpsbim$REF == susiesnpsbim$V6),]

    b <- susiesnpsbim[!(susiesnpsbim$REF=="A" & susiesnpsbim$ALT=="T") & !(susiesnpsbim$REF=="T" & susiesnpsbim$ALT=="A") & !(susiesnpsbim$REF=="C" & susiesnpsbim$ALT=="G") &
                        !(susiesnpsbim$REF=="G" & susiesnpsbim$ALT=="C") ]

    compmatchbim <- b[(b$A1.x==values(bases,keys=b$V5) & b$REF == values(bases,keys=b$V6)),]
    flippedbim <- b[(b$A1.x==b$V6 & b$REF == b$V5),]
    compflippedbim <- b[(b$A1.x==values(bases,keys=b$V6) & b$REF == values(bases,keys=b$V5)),]

    if(dim(flippedbim)[1] > 0){
      flippedbim <- mutate(flippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
    }
    if(dim(compflippedbim)[1] > 0){
      compflippedbim <- mutate(compflippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
    }

    matchbimsnps <- rbind(matchbim, compmatchbim, flippedbim, compflippedbim) %>% arrange(POS)

    snplist <- matchbimsnps$ID

    x <- genomat[,colnames(genomat) %in% snplist]
    R <- cor(x)

    susiesnps <- filter(matchbimsnps, ID %in% snplist)
    gwascolocsusie <- list("beta" = susiesnps$BETA.y, "varbeta" = (susiesnps$SE.y)^2, "snp" = susiesnps$ID, "position" = susiesnps$POS, "type" = "cc", "LD" = R, "N" = 100000)
    qtlcolocsusie <- list("beta" = setNames(susiesnps$BETA.x, susiesnps$ID), "varbeta" = (susiesnps$SE.x)^2, "snp" = susiesnps$ID, "position" = susiesnps$POS, "type" = "quant", "N" = susiesnps$OBS_CT[1], "MAF" = susiesnps$A1_FREQ, "sdY"=1, LD = R)
    saveRDS(gwascolocsusie, file = file.path(parent_dir, target_gene, paste0(target_gene, "_gwascolocsusie")))
    saveRDS(qtlcolocsusie, file = file.path(parent_dir, target_gene, paste0(target_gene, "_qtlcolocsusie")))
    saveRDS(R, file = file.path(parent_dir, target_gene, paste0(target_gene, "_LDcorr")))
  }
}