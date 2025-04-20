import os
import sys
import argparse

def check_arg(args=None):
    """
    Function to parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
    parser.add_argument("-r", "--reference",
    help="reference transcriptome",
    required=True)
    parser.add_argument("-c", "--condensed",
    help="condensed pipeline", 
    required=True)
    return parser.parse_args(args)



def checkLDdata(output_folder: str, ):
    print("MVA: Option specified to download LD data from the internet!")
    print("MVA: Obtaining LDdata (WARNING: This is a large download and will likely take time...)")
    if not os.path.isdir(f"{output_folder}/1000g_vcfs"):
        os.mkdir(f"{output_folder}/1000g_vcfs")
    
    os.chdir(output_folder)
    os.chdir("1000g_vcfs")
    
    print("MVA: Checking and downloading if any files are missing from LD directory")
    for x in range(1, 23):
        if not os.path.isfile(f"ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"):
            os.system(f"wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
        if not os.path.isfile(f"ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"):
            os.system(f"wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi")
        if not os.path.isfile(f"20130606_g1k.ped"):
            os.system(f"wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
    os.chdir("..")
    os.chdir("..")

def main():
    downloadLDdata("test_folder")

if __name__ == "__main__":
    main()