# CAP - coloc Automation Pipeline

# Introduction 
This pipeline provides an intuitive approach to colocalization, capable of processing ARIC-formatted pQTLs, GTEx-formatted eQTLs, and GWAS summary statistics. It tests for shared genetic signals between GWAS and QTL data by identifying overlapping causal variants, helping to pinpoint genes or features that are not only associated with a trait but also have their expression influenced by the same variants. Building on previous work, this pipeline allows users to input a list of genes and run colocalization analysis without modifying the script for each iteration.

## Dependencies 
* R Packages
  * data.table
  * dplyr
  * coloc
  * hash
  * optparse
  * locuscompareR
  * plink

## Description of Scripts 
* preprocess.R
  * Used for data cleanup and preprocessing. It parses command-line inputs, creates output directoreis, loads and prepares data into coloc format for downstream use.
* svassumption.R
  * Used to carry out single variant analysis.
* locuscompare.R
  * Used to carry out multiple variant analysis and visulation using LocusCompare.
* wrapper.py
  * Combines scripts.
 
## Specifications for Gene Input Data
* Using GTEx-formatted eQTLs
  * A .txt file with newline-seperated data. The first line is a header indicating the column content. Each subsequent line contains the gene using its Ensembl Gene ID (e.g. ENSG00000227232, ENSG00000162591, etc.).
* Using ARIC-formatted pQTLs
  * A .txt file with newline-seperated data. The first line is a header indicating the column content. Each subsequent line contains the gene using its gene symbol (e.g. LAYN, PTEN, TP53I3, etc.).
* See example files for reference (genes_eqtl.txt, genes_eqtl.txt)
 
## Necessary Inputs 
To run the script, the following input options are required. Ensure that each file and directory is correctly specified:
* `--process`: The type of data processing to be used. Default is "pqtl".
  * If you are using eQTL data, add the process flag `--process "eqtl"`, otherwise do not.
* `--genes`: Path to where list of genes of interest is located. For specific formatting instructions, refer to the "Specifications for Gene Input Data" section in the README.
* `--seqIDdir`: File path for the ARIC pQTLs seqid.txt file. Only needed if using pQTL data.
* `--pQTLdir`: Directory path where pQTL data is stored. Only needed if using pQTL data.
* `--eQTLdir`: Directory path where eQTL data is stored. Only needed if using eQTL data.
* `--GWASdir`: File path for the GWAS summary statistics data.
* `--CHR_input`: The column name in the GWAS summary statistics for chromosome number.
* `--BP_input`: The column name in the GWAS summary statistics for base pair position.
* `--A1_input`: The column name in the GWAS summary statistics for the effect allele (A1).
* `--A2_input`: The column name in the GWAS summary statistics for the alternate allele (A2).
* `--BETA_input`: The column name in the GWAS summary statistics for the BETA value.
* `--SE_input`: The column name in the GWAS summary statistics for the standard error (SE).
* `--outputdir`: Name of the directory where the output will be stored. The directory will be created automatically.

## Testing preprocess.R
* As the p/eQTL and GWAS summary statistic files are too large to upload here, please use the class server for testing as we figure out how to provide test data.
* Test Command for eQTL:
  * `Rscript /home/project4/preprocessV2.R --process "eqtl" --genes "/home/project4/genes_eqtl.txt" --GWASdir "/home/data/coloc_project/gwas_sumstats/breast_hg38.txt.gz" --eQTLdir "/home/project4/some_gtex" --CHR_input "CHR" --BP_input "BP" --A1_input "A1" --A2_input "A2" --BETA_input "BETA" --SE_input "SE" --outputdir "eQTL_results"`
* Test Command for pQTL:
  * `Rscript /home/project4/preprocessV2.R --genes "/home/project4/genes_pqtl.txt" --seqIDdir "/home/data/coloc_project/ARIC_pQTLs/seqid.txt" --GWASdir "/home/data/coloc_project/gwas_sumstats/breast_hg38.txt.gz" --pQTLdir "/home/data/coloc_project/ARIC_pQTLs/EA/" --CHR_input "CHR" --BP_input "BP" --A1_input "A1" --A2_input "A2" --BETA_input "BETA" --SE_input "SE" --outputdir "eQTL_results"`

## Testing locuscompare.R 
* `Rscript /home/project4/locuscompare.R --gwas ./eQTL_results/ENSG00000116254/ENSG00000116254_gwascoloc --twas ./eQTL_results/ENSG00000116254/ENSG00000116254_qtlcoloc --ref_dir /home/data/coloc_project/ref_data/ --gene_dir ./eQTL_results/ENSG00000116254/`
