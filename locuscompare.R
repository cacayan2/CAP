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
  make_option("--genes", type = "character", help = "List of genes used for the analysis"),
  make_option("--outputdir", type = "character", help = "Name of directory where you want output, this will be made for you")
)

opt <- parse_args(OptionParser(option_list = option_list))

genes <- read.table(opt$genes)
genes <- genes$V1
print(str(genes))
out_dir <- opt$outputdir

save.image("~/test.RData")
coloc_construct <- function(genes){
  for(gene_name in genes){
    print(str(gene_name))
    if(file.exists(paste0(out_dir, gene_name, "/", gene_name, "_gwascoloc")) & file.exists(paste0(out_dir, gene_name, "/", gene_name, "_qtlcoloc"))){
      gene.name <- gene_name
      gene.dir <- paste0(out_dir, gene_name, "/")
      system(paste0('touch ', gene.dir, gene_name, ".RData"))
      r.data <- paste0(gene.dir, gene_name, ".RData")
      if(!dir.exists(paste0(gene.dir, "Results"))){system(paste0('mkdir ', gene.dir, "Results"))}
      res.dir <- paste0(gene.dir, "Results/")
      gwascoloc <- readRDS(paste0(gene.dir, gene_name, "_gwascoloc"))
      qtlcoloc <- readRDS(paste0(gene.dir, gene_name, "_qtlcoloc"))
      sva <- readRDS(file = paste0(gene.dir, gene_name, "_sva"))
      system(paste0('touch ', res.dir, gene_name, "_sva.csv"))
      fwrite(sva$results, file = paste0(res.dir, gene_name, "_sva.csv"))
      gwascolocsusie <- readRDS(file = paste0(gene.dir, gene_name, "_gwascolocsusie"))
      if(is.null(check_dataset(gwascolocsusie,req="LD"))){
        message("Successfully loaded GWAS data with LD Matrix!")
      }else{
        stop("Could not load GWAS data with LD Matrix.")
      }
      qtlcolocsusie <- readRDS(file = paste0(gene.dir, gene_name, "_qtlcolocsusie"))
      if(is.null(check_dataset(qtlcolocsusie,req="LD"))){
        message("Successfully loaded QTL data with LD Matrix!")
      }else{
        stop("Could not load QTL data with LD Matrix.")
      }
      susiesnps <- readRDS(file = paste0(gene.dir, gene_name, "_susiesnps"))
      gwassnps.df <- dplyr::select(susiesnps, c('RS', 'P.x'))
      twassnps.df <- dplyr::select(susiesnps, c('RS', 'P.y'))
      colnames(gwassnps.df) <- c('rsid', 'pval'); colnames(twassnps.df) <- c('rsid', 'pval')
      gwassnps.id <- read_metal(gwassnps.df, marker_col = 'rsid', pval_col = 'pval')
      twassnps.id <- read_metal(twassnps.df, marker_col = 'rsid', pval_col = 'pval')
      merged <- merge(gwassnps.id, twassnps.id, by = 'rsid', suffixes = c('1','2'), all = FALSE)
      lead.snp <- get_lead_snp(merged)
      snp.pos <- get_position(merged, 'hg38')
      leadsnp.loc <- dplyr::filter(snp.pos,rsid == lead.snp)
      sgwas = runsusie(gwascolocsusie)
      if(!is.null(sgwas$sets$cs)){
        gwas_cc.count <- length(sgwas$sets$cs)
        gwas.count <- sum(sapply(sgwas$sets$cs, length))
        stwas = runsusie(qtlcolocsusie)
      }else{
        gwas.count <- 0
        gwas_cc.count <- 0
        twas.count <- 0
        twas_cc.count <- 0
        save(gene.name, gwascoloc, qtlcoloc, sva, gwascolocsusie, qtlcolocsusie, susiesnps, gwassnps.df, twassnps.df, lead.snp, sgwas,
             gwas_cc.count, gwas.count, twas_cc.count, twas.count, file = r.data)
        rmarkdown::render('./visualize.Rmd', params = list(rdata_path = r.data), output_file = paste0(res.dir, gene_name, "_coloc.html"))
        message("No credible sets detected for the GWAS data in this region for this phenotype; moving to the next gene.")
      }
      if(!is.null(stwas$sets$cs)){
        twas_cc.count <- length(stwas$sets$cs)
        twas.count <- sum(sapply(stwas$sets$cs, length))
        susie.res = coloc.susie(sgwas, stwas)
        susie.sum <- susie.res$summary
        colnames(susie.sum) <- c('nSNPs', 'GWAS SNP', 'TWAS SNP', 'H0', 'H1', 'H2','H3', 'H4', 'GCS', 'TCS')
        susie.sum <- as.data.frame(susie.sum)
        susie.full <- susie.res$results
        system(paste0('touch ', res.dir, gene_name, '_mva.csv'))
        fwrite(susie.full, file = paste0(res.dir, gene_name, '_mva.csv'))
        save(gene.name, gwascoloc, qtlcoloc, sva, gwascolocsusie, qtlcolocsusie, susiesnps, gwassnps.df, twassnps.df, lead.snp, sgwas,
             stwas, gwas_cc.count, gwas.count, twas_cc.count, twas.count,susie.res, susie.sum, file = r.data)
        rmarkdown::render(paste0('./visualize.Rmd'), params = list(rdata_path = r.data), output_file = paste0(res.dir, gene_name, "_coloc.html"))
      } else{
        twas.count <- 0
        twas_cc.count <- 0
        save(gene.name, gwascoloc, qtlcoloc, sva, gwascolocsusie, qtlcolocsusie, susiesnps, gwassnps.df, twassnps.df, lead.snp, sgwas,
             stwas, gwas_cc.count, gwas.count, twas_cc.count, twas.count, file = r.data)
        rmarkdown::render('./visualize.Rmd', params = list(rdata_path = r.data), output_file = paste0(res.dir, gene_name, "_coloc.html"))
        message("No credible sets detected for the QTL data in this region for this phenotype; moving to next gene.")
      }
    }else{
      message(paste0("The gene ", gene_name, " is not found in this region!"))
    }
    
  }
}

coloc_construct(genes)