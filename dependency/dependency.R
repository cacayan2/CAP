# The packages used in the CAP pipeline:
# [1] "devtools"           "locuscomparer"       "cgwtools"            "data.table"         
# [5] "dplyr"               "coloc"               "hash"                "optparse"           
# [9] "R.utils"             "stringr"             "parallel"

if(!require(devtools)){install.packages("devtools")}
if(!require(locuscomparer)){devtools::install_github("boxiangliu/locuscomparer")}
if(!require(cgwtools)){install.packages("cgwtools")}
if(!require(data.table)) {install.packages("data.table")}
if(!require(dplyr)) {install.packages("dplyr")}
if(!require(coloc)) {install.packages("coloc")}
if(!require(hash)) {install.packages("hash")}
if(!require(optparse)) {install.packages("optparse")}
if(!require(R.utils)) {install.packages("R.utils")}
if(!require(stringr)) {install.packages("stringr")}
if(!require(parallel)) {install.packages("parallel")}