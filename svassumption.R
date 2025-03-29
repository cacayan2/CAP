

# Loading appropriate libraries. 
invisible(library(data.table))
invisible(library(dplyr))
invisible(library(coloc))
invisible(library(hash))

# Obtain command line arguments.
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  print("Please pass an argument to the command.")
  stop()
}

superpop <- args[1]

# Steps 1-3, to be deleted!
"%&%" <- function(a, b) paste(a, b, sep = "")


pQTLdir <- "coloc_project/ARIC_pQTLs/EA/"
GWASdir <- "coloc_project/gwas_sumstats/"

protein <- "SeqId_2635_61"
cancer <- "breast"

print("Reading QTL and GWAS data...")
pQTLs = fread(pQTLdir %&% protein %&% ".PHENO1.glm.linear")
gwas = fread(GWASdir %&% cancer %&% "_hg38.txt.gz")

print("Joining QTL and gwas data...")
joined = inner_join(pQTLs, gwas, by = c("ID" = "RS"))

print("Adjusting for flipped and ambiguous bases...")
bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"

match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

a <- joined[!(joined$REF == "A" & joined$ALT == "T") & !(joined$REF == "T" & joined$ALT == "A") & !(joined$REF == "C" & joined$ALT == "G") & !(joined$REF == "G" & joined$ALT == "C")]

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]
flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

if(dim(flipped)[1] > 0){
  flipped <- mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped <- mutate(compflipped,BETA.y = -1*BETA.y)
}

matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

print("Generating coloc datasets...")
gwascoloc <- list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "cc")
pqtlcoloc <- list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )

print("Running coloc assuming single variant assumption...")
# Running coloc using single variant assumption
sv.res <- coloc.abf(dataset1 = gwascoloc, dataset2 = pqtlcoloc)
sv.sens <- sensitivity(sv.res, "H4 > 0.5")

print("Setting up linkage disequilibrium data...")
