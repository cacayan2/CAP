from imports import *

class LDcreate():
    """
    Class for creating LD data.
    User either provides their own LD data or uses this class to download LD data. 
    """

    @abstractmethod
    def check_args(args=None):
        """
        Function to parse command line arguments.
        Args:
            args (list, optional): List of arguments passed from the command line.
        Returns:
            argparse.Namespace: Parsed arguments object.
        """
        parser = argparse.ArgumentParser(description="ldcreate.py")
        parser.add_argument("-c", "--concurrent",
                            help="Number of conccurrent downloads for aria2c.", default="16",
                            required=False)
        parser.add_argument("-t", "--threads",
                            help="Number of threads to use.", 
                            required=True)
        return parser.parse_args(args)

    @abstractmethod
    def create_LD(args: argparse.Namespace):
        """
        Function to create LD data.
        """
        print("LD_SETUP: Creating directory for LD data...")
        if not os.path.dir("1000g_vcfs"):
            # Creating directories for LD data. 
            print("LD_SETUP: Creating directory for LD data!")
            os.mkdir("1000g_vcfs")
            os.chdir("1000g_vcfs")
            os.mkdir("raw")
            os.chdir("raw")
            
            # Starting download of required files. 
            print("LD_SETUP: Initiating download of LD data!")
            os.system(f"aria2c --allow-overwrite=true -x {args.concurrent} -m 0  -i urls.txt")
            os.chdir("..")
            os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")

            # Concatenating LD data.
            print("LD_SETUP: Concatenating LD data!")
            os.system(f"bcftools concat -0z --threads {args.threads} -o ALL.chr.vcf.gz raw/ALL.chr*.vcf.gz")

            # Indexing LD data.
            print("LD_SETUP: Indexing LD data!")
            os.system(f"tabix -p vcf ALL.chr.vcf.gz")

            # Removing raw LD data.
            print("LD_SETUP: Removing raw LD data!")
            os.system(f"rm -r raw")

            print("LD_SETUP: LD data creation complete!")

        else:
            SystemExit("LD_SETUP: Error, LD data folder already exists!")
        

def main():
    """
    Main function to handle the execution flow of the pipeline.
    Parses arguments, performs sanity checks, and runs preprocessing and LD data setup.
    """
    # Parse command line arguments.
    arguments = LDcreate.check_args(sys.argv[1:])

    # Create LD data.
    LDcreate.create_LD(arguments)

    