# CAP - `coloc` Automation Pipeline
## Authors: Emil Cacayan, Cian Dotson, Richa Patel
### Table of Contents
* **[Introduction](#introduction)**
* **[Data](#data)**
* **[CAP Workflow](#cap-workflow)**
* **[Usage](#usage)**
	* **[Dependencies](#dependencies)**
	* **[Input Specifications](#input-specifications)**
		* **[Specifications for Gene Input Data](#specifications-for-gene-input-data)**
		* **[Linkage Disequilibrium Data](#linkage-disequilibrium-data)**
	* **[The `config.ini` File](#the-config-ini-file)**
	* **[Command Arguments](#command-arguments)**
* **[References](#references)**

### Introduction
This pipeline provides an intuitive approach to colocalization, capable of processing ARIC-formatted pQTLs, GTEx-formatted eQTLs, and GWAS summary statistics. It tests for shared genetic signals between GWAS and QTL data by identifying overlapping causal variants, helping to pinpoint or features that are not only associated with a trait but also have their expression influenced by the same variants. Building on previous work, this pipeline allows users to input a list of genes and run colocalization analysis without modifying the script for each iteration.

### Data
This pipeline is a modification of a pipeline created to generate colocalization for pQTL data in a paper published detailing a proteome association study of breast, prostate, ovarian, and endometrial cancers (Gregga). The intended use of this pipeline is to be a general-use pipeline to run colocalization analyses across an array of data inputs for both eQTL and pQTl data, as opposed to the hard-coded inputs used in the study. The data provided for testing is a truncated version of the data used in the aforementioned study. 

### CAP Workflow

### Usage
Please note - all of the standard output for the scripts are sent to the `CAP.log` created in the directory upon running the script. To view progress of the pipeline or if any troubleshooting is necessary, please view the output in the log file. The log file is refreshed upon each run of the pipeline. 

#### Dependencies
This pipeline is designed for a Unix environment, and requires the following software to function:
* [Linux/Unix](https://www.linux.org/pages/download/)
* [Python3](https://www.python.org/downloads/)
* [plink2](https://www.cog-genomics.org/plink/2.0/)
* [bcftools](https://www.htslib.org/download/)

In addition to this software, the following `R` packages are utilized (which are automatically installed and loaded by the pipeline if not yet installed):
* [`data.table`](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
* [`dplyr`](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)
* [`coloc`](https://chr1swallace.github.io/coloc/)
* [`hash`](https://cran.r-project.org/web/packages/hash/index.html)
* [`optparse`](https://cran.r-project.org/web/packages/optparse/index.html)
* [`R.utils`](https://cran.r-project.org/web/packages/R.utils/index.html)
* [`stringr`](https://cran.r-project.org/web/packages/stringr/index.html)

#### Input Specifications
##### Specifications for Gene Input Data
As mentioned, this pipeline accepts ARIC-formatted pQTLs or GTEx-formatted eQTLs and GWAS summary statistics to perform colocalization. Due to the lack of a consensus in the format of GWAS, eQTL, and pQTL data, we ask that the following information be included in GWAS, eQTL, and pQTL data with well-defined column names, entered in the [`config.ini`](#the-config-ini-file) file:

* Chromosome number
* Base pair position
* Effect Allele
* Alternate allele
* BETA value
* Standard Error 

This pipeline runs two colocalization procedures with the `coloc` package in `R`, the first being the assumption of 0 or 1 causal variant in each trait (single variant assumption), and the other is the understanding that multiple causal variants can be involved in the shared genetic influence between two traits (multiple variant assumption). The former is a trivial task, but the latter requires linkage disequilibrium data (represented as a correlation matrix between SNP's) to help cluster tightly linked variants and decrease artificially inflated false positive rates. 

The linkage disequilibrium data is generated from data obtained from phase 3 of the 1000 genomes project, a project aimed to map the majority of human genetic variation. The data that is required for this pipeline is extremely large, and we highly recommend that upon downloading the user keeps the data in a safe place for reuse. This pipeline has the added capability of downloading the required data built-in. The pipeline can also accommodate already downloaded and processed data skipping corresponding steps, streamlining the process for future iterations of the pipeline's usage. 

The output of this pipeline is a directory with a text file containing the associated $p$-values for the colocalization, and visualization of different plots via `RMarkdown` utilizing `locuscompareR` and other visualization tools. This tool is primarily intended to generate raw data, and does not have robust capability to interpret results. 

In order to select regions of interest, the pipeline requires a list of genes to generate an SNP list from for either eQTL or pQTL data. The pipeline requires the user specify a directory containing a list of such identifiers/gene symbols.
* **Using GTEx-formatted eQTLs**
	* A `.txt` file with newline-separated data. The first line is a header indicating the column content. Each subsequent line contains the gene using its Ensembl Gene ID (e.g. ENSG00000227232, ENSG00000162591, etc.).
* **Using ARIC-formatted pQTLs**
	* A `.txt` file with newline-separated data. The first line is a header indicating the column content. Each subsequent line contains the gene using its gene symbol (e.g. LAYN, PTEN, TP53I3, etc.).

##### The `config.ini` File
To reduce the verboseness of running the pipeline, a `config.ini` file is included in the repository. This file stores user-specified information for directories and preferences for the pipeline. Each section has the following keys which must be assigned values as followed:
* **`genes`**: Path to where the list of genes of interest is located. For specific formatting instructions, refer to the [Specifications for Gene Input Data](#specifications-for-gene-input-data) section.
* **`seqIDdir`**: File path for the ARIC pQTLs `seqid.txt` file. Only needed if using pQTL data.
* **`data_dir`**: Directory where either eQTL or pQTL data is stored. 
* **`CHR_input`**: The column name in the GWAS summary statistics for chromosome number.
* **`BP_input`**: The column name in the GWAS summary statistics for base pair position.
* **`A1_input`**: The column name in the GWAS summary statistics for the effect allele (A1).
* **`A2_input`**: The column name in the GWAS summary statistics for the alternate allele (A2).
* **`BETA_input`**: The column name in the GWAS summary statistics for the regression coefficient (BETA).
* **`SE_input`**: The column name in the GWAS summary statistics for standard error (SE). 
* **`LD_dir`**: Location of the LD data directory if not using pipeline to download from the internet. This will be ignored if the command line flag `-l`/`--lddownload` is set to true. 

When specifying a directory (not a file), please denote with an extra `/` at the very end of the path, for example (a directory called 1000g_vcfs):

```
LD_dir = "data/1000g_vcfs/"
```

##### Linkage Disequilibrium Data
The longest step involved in this pipeline is obtaining and processing the linkage disequilibrium data required for colocalization under the multiple variant hypothesis. There are a number of ways the pipeline streamlines this process.

To create the LD correlation matrix, the corresponding `.vcf` files for all chromosomes must be downloaded from the [1000 Genomes Project phase 3 release](https://www.internationalgenome.org/data-portal/data-collection/phase-3). This pipeline has the capability to do this but the user can also download the data directly from [this directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) using the following command:

```
for chr in {1..22}; do
  wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt
```

Save these to a directory and specify the directory in the [`config.ini`](#the-config-ini-file) under the `LD_dir` section. When running the `wrapper.py` script, please specify the flag `-l/--lddownload` as false, otherwise the directory in the `config.ini` file will be disregarded. The pipeline will verify that all files are present in the directory (not if they have been corrupted or are incomplete). 

Otherwise, to make the pipeline download the data, simply run with the `-l/--lddownload` flag as true,  and a directory containing the data will be created for you in a directory called `1000g_vcfs`. If the folder is created and is complete, the pipeline will not download the data again and skip to processing.

The processing step is also non-trivial, and the user can specify the number of threads to run processing with using the `-t/--threads` flag. This step involves the creation of a new directory called `LD_output`. Again, if this folder is present and complete (i.e. a round of processing on the downloaded data is already finished), then this step will not be repeated. The output of all processing is a correlation matrix to be used for colocalization.

##### Command Arguments
To run the script, please clone the repository:
```
git clone https://github.com/cacayan2/CAP.git
```
The repository also contains a test data set that can be used to troubleshoot the pipeline. All of the scripts are managed by `wrapper.py`.

To run the pipeline, use the following command:

```
python3 wrapper.py -p <eqtl/pqtl> -t <t/f> -o <output> -s <superpop> -l <ld> -c <threads>
```
Here is a brief explanation of all the command line arguments:

* `-p/--process`: Specify whether using eQTL or pQTL data.
	* Acceptable values: "eqtl", "pqtl"
* `-t/--test`: Specify whether using test data or real data.
	* Acceptable values: "True", "true", "T", "t", or "False", "false", "F", "f"
* `-s/--superpop`: The subset of the population to use for generation of the linkage disequilibrium data. 
	* Acceptable values: Please see below for a list of available populations.
* `-l/--ld`: Specifies whether the pipeline will download the linkage disequilibrium data or if the user is providing linkage disequilibrium data. 
	* Acceptable values: "True", "true", "T", "t", or "False", "false", "F", "f"
* `-c/--threads`: The number of threads used for processing - ideal number of threads is environment dependent.
	* Acceptable values: Whole number values up to the number of cores available on workspace.

### References
1. Gregga I, Pharoah PDP, Gayther SA, Manichaikul A, Im HK, Kar SP, Schildkraut JM, Wheeler HE. Predicted Proteome Association Studies of Breast, Prostate, Ovarian, and Endometrial Cancers Implicate Plasma Protein Regulation in Cancer Susceptibility. Cancer Epidemiol Biomarkers Prev. 2023 Sep 1;32(9):1198-1207. doi: 10.1158/1055-9965.EPI-23-0309. PMID: 37409955; PMCID: PMC10528410.
