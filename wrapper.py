import os
import sys
import argparse
import configparser
from constants import *
from distutils.util import strtobool

def check_arg(args=None):
    """
    Function to parse command line arguments.
    Args:
        args (list, optional): List of arguments passed from the command line.
    Returns:
        argparse.Namespace: Parsed arguments object.
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
    parser.add_argument("-s", "--superpop", default="ALL",
                        help="Name of superpopulation to pull LD data from", required=True)
    parser.add_argument("-l", "--ld", default="True",
                        help="Download LD data from the internet", required=True)
    parser.add_argument("-c", "--threads", default="2",
                        help="Number of threads to use", required=True)
    
    return parser.parse_args(args)

def sanity_check(process, test, superpop, ld):
    """
    Perform sanity checks on command line arguments to ensure they are valid.
    Args:
        process (str): The type of data processing ('eqtl' or 'pqtl').
        test (str): The test flag ('True' or 'False').
        superpop (str): The superpopulation flag (e.g., 'ALL', 'EUR', 'AFR', etc.).
        ld (str): The LD flag ('True' or 'False').
    """
    # Check for valid process type
    if process != "eqtl" and process != "pqtl":
        raise SystemExit("PREPROCESS: Process flag must be either 'eqtl' or 'pqtl'")
    
    # Check for valid test flag
    if test not in ["True", "true", "False", "false", "T", "t", "F", "f"]:
        raise SystemExit("PREPROCESS: Test flag must be True or False.")
    
    # Check for valid superpopulation flag
    if superpop not in ["ALL", "EUR", "AFR", "AMR", "EAS", "SAS", 
                        "YRI", "LWK", "GWD", "MSL", "ESN", "ASW",
                        "ACB", "MXL", "PUR", "CLM", "PEL", "GIH", 
                        "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", 
                        "TSI", "FIN", "IBS", "GBR", "IBS", "GIH",
                        "PJL", "BEB", "STU", "ITU"]:
        raise SystemExit("PREPROCESS: Superpop flag must be one of 'ALL', 'EUR', 'AFR', 'AMR', 'EAS', or 'SAS'")
    
    # Check for valid LD flag
    if ld not in ["True", "true", "False", "false", "T", "t", "F", "f"]:
        raise SystemExit("PREPROCESS: LD flag must be True or False.")

def checkLDdata(output_folder: str, ld: str):
    """
    Check and set up LD data directory.
    Args:
        output_folder (str): The folder where LD data should be stored.
        ld (str): Flag indicating whether to download new LD data.
    Returns:
        str: Path to the LD data directory.
    """
    current_wd = os.getcwd()
    ld_dir = ""
    
    # If LD data needs to be downloaded
    if ld in ["True", "true", "T", "t"]:
        print("LD_SETUP: Creating directory for LD data...")
        if not os.path.isdir("1000g_vcfs"):
            os.mkdir("1000g_vcfs")
        os.chdir("1000g_vcfs")
        ld_dir = "1000g_vcfs/"
    else:
        # If LD data folder already exists, use it
        if not os.path.isdir(output_folder):
            SystemExit("LD_SETUP: LD data folder specified does not exist. Please create folder and try again.")
        ld_dir = output_folder
        os.chdir(ld_dir.replace("/", "").replace("\"", ""))

    print("LD_SETUP: Downloading LD data if needed...")
    # Download required files if they don't exist
    for x in range(1, 23):
        if not os.path.isfile(f"ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"):
            os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{str(x)}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    if not os.path.isfile(f"integrated_call_samples_v3.20200731.ALL.ped"):
        os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")
    if not os.path.isfile(f"GCF_000001405.40.gz"):
        os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress -c https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz")
    if not os.path.isfile("GCF_000001405.40.gz.tbi"):
        os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress -c https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi")

    os.chdir(current_wd)
    return ld_dir

def main():
    """
    Main function to handle the execution flow of the pipeline.
    Parses arguments, performs sanity checks, and runs preprocessing and LD data setup.
    """
    # Parse command-line arguments
    arguments = check_arg(sys.argv[1:])
    
    # Checking if arguments are correct. 
    print("WRAPPER: Verifying arguments are correct...")
    sanity_check(arguments.process, arguments.test, arguments.superpop, arguments.ld)
    
    # Determine if using test data or real data
    if arguments.test in ["True", "true", "T", "t"]:
        print("WRAPPER: Using test data!")
    else:
        print("WRAPPER: Using real data!")

    constants = Constants(arguments.test)

    print("WRAPPER: Creating log file...")
    # Remove existing log file and create a new one
    if os.path.exists("CAP.log"):
        os.system("rm CAP.log")
    with open("CAP.log", "w") as f:
        f.write("---Colocalization Automation Pipeline Log File---")
    
    # Preprocess data based on the process type (eQTL or pQTL)
    if arguments.process == "eqtl":
        print("WRAPPER: Preprocessing eQTL data...")
        os.system(f"Rscript preprocess.R --process {arguments.process} --genes {constants.genes()} --GWASdir {constants.GWASdir()} --eQTLdir {constants.data_dir()} --CHR_input {constants.CHR_input()} --BP_input {constants.BP_input()} --A1_input {constants.A1_input()} --A2_input {constants.A2_input()} --BETA_input {constants.BETA_input()} --SE_input {constants.SE_input()} --outputdir {arguments.output} >> CAP.log 2>&1")
    else:
        print("WRAPPER: Preprocessing pQTL data...")
        os.system(f"Rscript preprocess.R --process {arguments.process} --genes {constants.genes()} --seqIDdir {constants.seqID_dir()} --GWASdir {constants.GWASdir()} --pQTLdir {constants.data_dir()} --CHR_input {constants.CHR_input()} --BP_input {constants.BP_input()} --A1_input {constants.A1_input()} --A2_input {constants.A2_input()} --BETA_input {constants.BETA_input()} --SE_input {constants.SE_input()} --outputdir {arguments.output} >> CAP.log 2>&1")
    
    print("WRAPPER: Completed preprocessing! Generating LD data for multivariant assumption colocalization...")
    
    # Handle LD data downloading and processing
    if arguments.ld in ["True", "true", "T", "t"]:
        ld_dir = checkLDdata("1000g_vcfs/", arguments.ld)
    else:
        ld_dir = checkLDdata(constants.LD_dir(), arguments.ld)
    
    # Run single variant assumption and LD data processing
    print("WRAPPER: Downloaded LD data! Starting single variant assumption analysis and LD data processing...")
    os.system(f"Rscript svassumption.R --process {arguments.process} --superpop {arguments.superpop} --output {arguments.output} --data {arguments.output} --lddir {ld_dir} --lddownload {arguments.ld} --threads {arguments.threads} >> CAP.log 2>&1")

    # Run locuscompare and visualization
    print("WRAPPER: Finished generating LD data processing and single variant assumption! Starting locuscompare...")
    os.system(f"Rscript locuscompare.R --genes {constants.genes()} --outputdir {arguments.output}/ >> CAP.log 2>&1")

if __name__ == "__main__":
    main()
    print("WRAPPER: Pipeline complete! See CAP.log for details.")