if(!require(devtools)){install.packages("devtools")}
suppressPackageStartupMessages(invisible(library(devtools)))
if(!require(locuscomparer)){devtools::install_github("boxiangliu/locuscomparer")}
if(!require(cgwtools)){install.packages("cgwtools")}


suppressPackageStartupMessages(invisible(library(data.table)))
suppressPackageStartupMessages(invisible(library(dplyr)))
suppressPackageStartupMessages(invisible(library(coloc)))
suppressPackageStartupMessages(invisible(library(hash)))
suppressPackageStartupMessages(invisible(library(optparse)))
suppressPackageStartupMessages(invisible(library(locuscomparer)))
suppressPackageStartupMessages(invisible(library(stringr)))
suppressPackageStartupMessages(invisible(library(cgwtools)))


# takes the binary files for the gwas and e/pqtl coloc objects as inputs for the Rscript #
option_list <- list(
  make_option("--genes", type = "character", help = "List of genes used for the analysis")
)

opt <- parse_args(OptionParser(option_list = option_list))

genes <- opt$genes 


coloc_construct <- function(genes){
  for(gene_name in genes){
    if(file.exists(paste0(getwd(), "/", gene_name, "/", gene_name, "_gwascoloc")) & file.exists(paste0(getwd(), "/", gene_name, "/", gene_name, "_qtlcoloc"))){
      gene.dir <- paste0(getwd(), "/", gene_name, "/")
      system(paste0('touch ', gene.dir, gene_name, ".RData"))
      r.data <- paste0(gene.dir, gene_name, ".RData")
      gwascoloc <- readRDS(paste0(gene.dir, gene_name, "_gwascoloc"))
      qtlcoloc <- readRDS(paste0(gene.dir, gene_name, "_qtlcoloc"))
      save(gwascoloc, qtlcoloc, file = r.data)
      gwascolocsusie <- readRDS(paste0(gene.dir, gene_name, "_gwascolocsusie"))
      if(is.null(check_dataset(gwascolocsusie,req="LD"))){
        message("Successfully loaded GWAS data with LD Matrix!")
      }else{
          stop("Could not load GWAS data with LD Matrix.")
        }
      resave(gwascolocsusie, file = r.data)
      qtlcolocsusie <- readRDS(paste0(gene.dir, gene_name, "_qtlcolocsusie"))
      if(is.null(check_dataset(qtlcolocsusie,req="LD"))){
        message("Successfully loaded QTL data with LD Matrix!")
      }else{
        stop("Could not load QTL data with LD Matrix.")
      }
      resave(qtlcolocsusie, file = r.data)
      susiesnps <- readRDS(paste0(gene.dir, gene_name, "_susiesnps"))
      resave(susiesnps, file = r.data)
      gwassnps.df <- dplyr::select(susiesnps, c('RS', 'P.x'))
      twassnps.df <- dplyr::select(susiesnps, c('RS', 'P.y'))
      colnames(gwassnps.df) <- c('rsid', 'pval'); colnames(twassnps.df) <- c('rsid', 'pval')
      resave(gwassnps.df, file = r.data)
      resave(twassnps.df, file = r.data)
      gwassnps.id <- read_metal(gwassnps.df, marker_col = 'rsid', pval_col = 'pval')
      twassnps.id <- read_metal(twassnps.df, marker_col = 'rsid', pval_col = 'pval')
      merged <- merge(gwassnps.id, twassnps.id, by = 'rsid', suffixes = c('1','2'), all = FALSE)
      lead.snp <- get_lead_snp(merged)
      resave(lead.snp, file = r.data)
      leadsnp.loc <- dplyr::filter(snp.pos,rsid == lead.snp)
      resave(leadsnp.loc, file = r.data)
      snp.pos <- get_position(merged, 'hg38')
      sgwas = runsusie(gwassusiecoloc)
      resave(sgwas, file = r.data)
      if(!is.null(sgwas$sets$cs)){
        gwas_cc.count <- length(sgwas$sets$cs)
        resave(gwas_cc.count, file = r.data)
        gwas.count <- sum(sapply(sgwas$sets$cs, length))
        resave(gwas.count, file = r.data)
        stwas = runsusie(twassusiecoloc)
        resave(stwas, file = r.data)
      }else{
        gwas.count <- 0
        resave(gwas.count, file = r.data)
        gwas_cc.count <- 0
        resave(gwas_cc.count, file = r.data)
        twas.count <- 0
        resave(twas.count, file = r.data)
        twas_cc.count <- 0
        resave(twas_cc.count, file = r.data)
        rmarkdown::render('./visualize.Rmd', params = list(rdata_path = r.data), output_file = paste0(gene_name, "_coloc.html"))
      }
      if(!is.null(stwas$sets$cs)){
        twas_cc.count <- length(stwas$sets$cs)
        resave(twas_cc.count, file = r.data)
        twas.count <- sum(sapply(stwas$sets$cs, length))
        resave(twas.count, file = r.data)
        susie.res = coloc.susie(sgwas, stwas)
        resave(susie.res, file = r.data)
        susie.sum <- susie.res$summary
        colnames(susie.sum) <- c('nSNPs', 'GWAS SNP', 'TWAS SNP', 'H0', 'H1', 'H2','H3', 'H4', 'GCS', 'TCS')
        susie.sum <- as.data.frame(susie.sum)
        resave(susie.sum, file = r.data)
        susie.full <- susie.res$results
        fwrite(susie.full, file = paste0(gene.dir, gene_name, '_susie.txt'))
        rmarkdown::render(paste0('./visualize.Rmd'), params = list(rdata_path = r.data), output_file = paste0(gene.dir, gene_name, "_coloc.html"))
      } else{
        twas.count <- 0
        resave(twas.count, file = r.data)
        twas_cc.count <- 0
        resave(twas_cc.count, file = r.data)
        rmarkdown::render('./visualize.Rmd', params = list(rdata_path = r.data), output_file = paste0(gene.dir, gene_name, "_coloc.html"))
      }
    }else{
      message(paste0("The gene ", gene_name, " is not found in this region!"))
    }
    
}
}