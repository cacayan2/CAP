# Testing preprocess.R (Steps 1-3)

Test Command for eQTL:

`Rscript /home/project4/DataPrepv2.R --process "eqtl" --genes "/home/project4/genes_eqtl.txt" --GWASdir "/home/data/coloc_project/gwas_sumstats/breast_hg38.txt.gz" --eQTLdir "/home/project4/some_gtex" --CHR_input "CHR" --BP_input "BP" --A1_input "A1" --A2_input "A2" --BETA_input "BETA" --SE_input "SE"`

Test Command for pQTL:

`Rscript /home/project4/DataPrepv2.R --genes "/home/project4/genes_pqtl.txt" --seqIDdir "/home/data/coloc_project/ARIC_pQTLs/seqid.txt" --GWASdir "/home/data/coloc_project/gwas_sumstats/breast_hg38.txt.gz" --pQTLdir "/home/data/coloc_project/ARIC_pQTLs/EA/" --CHR_input "CHR" --BP_input "BP" --A1_input "A1" --A2_input "A2" --BETA_input "BETA" --SE_input "SE"`


