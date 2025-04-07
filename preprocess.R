library(data.table)
library(dplyr)
library(coloc)
library(hash)
library(optparse)
library(R.utils)

# Define string concatenation function
"%&%" <- function(a, b) paste(a, b, sep="")

# Define command-line options
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--genes", type = "character", help = "List of genes of intrest, see READ.ME for specifications."),
  make_option("--seqIDdir", type = "character", help = "File path for ARIC pQTLs seqid.txt"),
  make_option("--GWASdir", type = "character", help = "File path for GWAS summary statistics"),
  make_option("--pQTLdir", type = "character", help = "Directory path for pQTL data"),
  make_option("--eQTLdir", type = "character", help = "Directory path for eQTL data"),
  make_option("--CHR_input", type = "character", help = "Column name for GWAS chromosome number"),
  make_option("--BP_input", type = "character", help = "Column name for GWAS bp postion"),
  make_option("--A1_input", type = "character", help = "Column name for GWAS affect allele (A1)"),
  make_option("--A2_input", type = "character", help = "Column name for GWAS alternate allele (A2)"),
  make_option("--BETA_input", type = "character", help = "Column name for GWAS BETA"),
  make_option("--SE_input", type = "character", help = "Column name for GWAS SE"),
  make_option("--outputdir", type = "character", help = "Name of directory where you want output, this will be made for you")
)

opt <- parse_args(OptionParser(option_list=option_list))

#gets working directory and the new dir the user wants to make for output
parent_dir <- file.path(getwd(), opt$outputdir) 
dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)  #create if it doesn't exist

#read in the genes of intrest
genes = fread(opt$genes)

#loop through each gene and create a directory
for (gene in genes$genes) {
  dir_path <- file.path(parent_dir, gene)
  dir.create(dir_path) 
}

#if using eqtl data
if (opt$process == "eqtl") {
  #read in gwas and eqtl datasets
  gwas = fread(opt$GWASdir)
  eQTLs <- fread(opt$eQTLdir)
  
  #format eqtl data, break up variant_id column and rename certain columns for better handling
  eQTLs <- eQTLs %>%
    mutate(`#CHROM` = sapply(strsplit(variant_id,"_"), `[`, 1), #chromosome number
           POS = sapply(strsplit(variant_id,"_"), `[`, 2), #position
           REF = sapply(strsplit(variant_id,"_"), `[`, 3), #reference allele, w/ indels
           A1 = sapply(strsplit(variant_id,"_"), `[`, 4)) %>% #alternate allele w/ indels
    rename(OBS_CT = ma_count, #number of allele observations
           A1_FREQ = maf, #A1 allele frequency
           BETA = slope, #regression coefficient for A1 allele
           SE = slope_se, #standard error
           P = pval_nominal) #p-value
  
  #filter eqtl data
  eQTLs <- eQTLs %>%
    filter(nchar(A1)==1, #remove indels from alt alleles
           nchar(REF) == 1) %>% #remove indels from ref alleles
    mutate(CHR = gsub("chr", "", `#CHROM`), #create chr column
           ALT = A1) #create ALT column
  eQTLs <- eQTLs %>%
    mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
  
  #remame GWAS input for better handling
  gwas <- gwas %>%
    rename_with(~ c("CHR", "POS", "A1", "A2", "BETA", "SE"), 
                .cols = c(opt$CHR_input, opt$BP_input, opt$A1_input, opt$A2_input, opt$BETA_input, opt$SE_input))
  gwas <- gwas %>%
    mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
  
  #for each gene in the target genes list, perform the following:
  for (target_gene in genes$genes) {
    gene_eQTL <- eQTLs %>% filter(grepl(paste0("^", target_gene), gene_id)) #filter eqtl data for specific gene
    if (nrow(gene_eQTL) == 0) {
      cat("Skipping gene", target_gene, "as it is not found in eQTL data\n") #if gene is not in data, skip it
      next
    }
    
    joined = inner_join(gene_eQTL,gwas, by=c("chr_pos"="chr_pos")) #join gwas and eqtl data based on chr_pos
    
    #check for complement bases
    #build hash table (like a Python dictionary)
    bases = hash()
    bases[["A"]] <- "T"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"
    bases[["T"]] <- "A"
    
    #pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
    match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]
    
    #remove ambiguous strand SNPs (A/T or C/G)
    a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
                 !(joined$REF=="G" & joined$ALT=="C") ]
    
    #of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)
    compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]
    
    flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
    compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]
    
    #if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
    if(dim(flipped)[1] > 0){
      flipped = mutate(flipped,BETA.y = -1*BETA.y)
    }
    if(dim(compflipped)[1] > 0){
      compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
    }
    
    #bind all and sort by position
    matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(chr_pos)
    
    matchsnps <- matchsnps %>%
      distinct(POS.y, .keep_all = TRUE) #keep only distinct
    matchsnps <- matchsnps %>%
      na.omit() #omit NAs
    
    #format gwas data for coloc analysis
    gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.y,
                     "type" = "cc")
    #format eqtl data for coloc analysis
    eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.x,
                     "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )
  
    #save the gwas/eqtl coloc formatted data in a binary file for downstream analysis
    saveRDS(gwascoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_gwascoloc")))
    saveRDS(eqtlcoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_qtlcoloc")))
    
    #save matchsnp data for LD downstream   
    saveRDS(matchsnps, file = file.path(parent_dir, target_gene, paste0(target_gene, "_matchsnps")))
    
  }
  #if --eqtl is not specified, default is pqtl 
} else {
  #read in necessary files
  gwas = fread(opt$GWASdir)
  pQTLdir = opt$pQTLdir
  seqid = fread(opt$seqIDdir)
  
  #remame GWAS input for better handling
  gwas <- gwas %>%
    rename_with(~ c("CHR", "POS", "A1", "A2", "BETA", "SE"), 
                .cols = c(opt$CHR_input, opt$BP_input, opt$A1_input, opt$A2_input, opt$BETA_input, opt$SE_input))
  gwas <- gwas %>%
    mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
  
  for (target_gene in genes$genes) {
    # Find the associated seqid_in_sample
    protein <- seqid %>%
      filter(entrezgenesymbol == target_gene) %>%
      pull(seqid_in_sample) #find the seqID which will be used to pull pqtl data
    
    if (length(protein) == 0) {
      cat("Skipping gene", target_gene, "as it is not found in seqID mapping\n")
      next
    }
    
    pQTLs = fread(pQTLdir %&% protein %&% '.PHENO1.glm.linear') #read in pqtl data for protein 
    
    #rename/format pqtl data
    pQTLs <- pQTLs %>%
      rename(CHR = `#CHROM`,
             POS = POS,
             REF = REF,
             ALT = ALT,
             A1 = A1, 
             BETA = BETA,
             SE = SE,
             OBS_CT = OBS_CT, 
             A1_FREQ = A1_FREQ)
    pQTLs <- pQTLs %>%
      mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
    
    joined = inner_join(pQTLs,gwas, by=c("chr_pos"="chr_pos")) #join gwas/pqtl data
    
    #check for complement bases
    #build hash table (like a Python dictionary)
    bases = hash()
    bases[["A"]] <- "T"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"
    bases[["T"]] <- "A"
    
    #pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
    match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]
    
    #remove ambiguous strand SNPs (A/T or C/G)
    a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
                 !(joined$REF=="G" & joined$ALT=="C") ]
    
    #of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)
    compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]
    flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
    compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]
    
    #if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
    if(dim(flipped)[1] > 0){
      flipped = mutate(flipped,BETA.y = -1*BETA.y)
    }
    if(dim(compflipped)[1] > 0){
      compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
    }
    
    #bind all and sort by position
    matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(chr_pos)
    
    #format gwas data for coloc analysis
    gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.y,
                     "type" = "cc")
    #format pqtl data for coloc analysis
    pqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.x,
                     "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )
    #save the gwas/eqtl coloc formatted data in a binary file for downstream analysis
    saveRDS(gwascoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_gwascoloc")))
    saveRDS(pqtlcoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_qtlcoloc")))
    
    #save matchsnp data for LD downstream 
    saveRDS(matchsnps, file = file.path(parent_dir, target_gene, paste0(target_gene, "_matchsnps")))
    
    
  }
  
}
