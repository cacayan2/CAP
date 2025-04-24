import os
import sys
import argparse
import configparser
from constants import *
from distutils.util import strtobool

def check_arg(args=None):
    """
    Function to parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
    parser.add_argument("-p", "--process", default="pqtl",
    help="If using eqtl data, add process flag, otherwise pqtl is the default.",
    required=True)
    parser.add_argument("-t", "--test", default="False",
    help="Test the pipeline with test data.", 
    required=True)
    parser.add_argument("-o", "--output",
    help="Name of directory where you want output, this will be made for you",
    required=True)
    parser.add_argument("-s", "--superpop",default="ALL",
    help="Name of superpopulation to pull LD data from", required=True)
    parser.add_argument("-l", "--ld",default="True",
    help="Download LD data from the internet", required=True)
    parser.add_argument("-c", "--threads", default="2",
    help="Number of threads to use", required=True)
    return parser.parse_args(args)

def sanity_check(process, test, superpop, ld):
    if process != "eqtl" and process != "pqtl":
        raise SystemExit("PREPROCESS: Process flag must be either 'eqtl' or 'pqtl'")
    if test != "True" and test != "true" and test != "False" and test != "false" and test != "T" and test != "t" and test != "F" and test != "f":
        raise SystemExit("PREPROCESS: Test flag must be True or False.")
    if superpop != "ALL" and superpop != "EUR" and superpop != "AFR" and superpop != "AMR" and superpop != "EAS" and superpop != "SAS":
        raise SystemExit("PREPROCESS: Superpop flag must be one of 'ALL', 'EUR', 'AFR', 'AMR', 'EAS', or 'SAS'")
    if ld != "True" and ld != "true" and ld != "False" and ld != "false" and ld != "T" and ld != "t" and ld != "F" and ld != "f":
        raise SystemExit("PREPROCESS: LD flag must be True or False.")

def checkLDdata(output_folder: str, ld: str):
    current_wd = os.getcwd()
    ld_dir = ""
    if ld == "True" or ld == "true" or ld == "T" or ld == "t":
        print("LD_SETUP: Creating directory for LD data...")
        if not os.path.isdir("1000g_vcfs"):
            os.mkdir("1000g_vcfs")
            os.chdir("1000g_vcfs")
        ld_dir = "1000g_vcfs/"
    else:
        ld_dir = output_folder

    os.chdir(ld_dir.replace("/", "").replace("\"", ""))

    print("LD_SETUP: Downloading LD data if needed...")
    for x in range(1, 23):
        if not os.path.isfile(f"ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"):
            os.system(f"wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    if not os.path.isfile(f"20130606_sample_info.txt"):
            os.system(f"wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt")

    os.chdir(current_wd)
    return ld_dir
def main():
    arguments = check_arg(sys.argv[1:])
    
    print("WRAPPER: Verifying arguments are correct...")
    sanity_check(arguments.process, arguments.test, arguments.superpop, arguments.ld)
    
    if not arguments.test == "True" or arguments.test == "true" or arguments.test == "T" or arguments.test == "t":
        print("WRAPPER: Using real data!")
    else:
        print("WRAPPER: Using test data!")

    constants = Constants(arguments.test)

    print("WRAPPER: Creating log file...")
    if os.path.exists("CAP.log"):
        os.system("rm CAP.log")
    with open("CAP.log", "w") as f:
        f.write("---Colocalization Automation Pipeline Log File---")
    
    if arguments.process == "eqtl":
        print("WRAPPER: Preprocessing eQTL data...")
        # os.system(f"Rscript preprocess.R --process {arguments.process} --genes {constants.genes()} --GWASdir {constants.GWASdir()} --eQTLdir {constants.data_dir()} --CHR_input {constants.CHR_input()} --BP_input {constants.BP_input()} --A1_input {constants.A1_input()} --A2_input {constants.A2_input()} --BETA_input {constants.BETA_input()} --SE_input {constants.SE_input()} --outputdir {arguments.output} >> CAP.log 2>&1")
    else:
        print("WRAPPER: Preprocessing pQTL data...")
        # os.system(f"Rscript preprocess.R --process {arguments.process} --genes {constants.genes()} --seqIDdir {constants.seqID_dir()} --GWASdir {constants.GWASdir()} --pQTLdir {constants.data_dir()} --CHR_input {constants.CHR_input()} --BP_input {constants.BP_input()} --A1_input {constants.A1_input()} --A2_input {constants.A2_input()} --BETA_input {constants.BETA_input()} --SE_input {constants.SE_input()} --outputdir {arguments.output} >> CAP.log 2>&1")
    
    print("WRAPPER: Completed preprocessing! Generating LD data for multivariant assumption colocalization...")
    
    if(arguments.ld == "True" or arguments.ld == "true" or arguments.ld == "T" or arguments.ld == "t"):
        ld_dir = checkLDdata("1000g_vcfs/", arguments.ld).replace("/", "")
    else:
        ld_dir = checkLDdata(constants.LD_dir(), arguments.ld)
    
    print("WRAPPER: Downloaded LD data! Starting single variant assumption analysis and LD data processing...")
    # os.system(f"Rscript svassumption.R --process {arguments.process} --superpop {arguments.superpop} --output {arguments.output} --lddir {ld_dir} --lddownload {arguments.ld} --threads {arguments.threads} >> CAP.log 2>&1")

if __name__ == "__main__":
    main()