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
        if not os.path.isdir("1000g_vcfs"):
            # Creating directories for LD data. 
            print("LD_SETUP: Creating directory for LD data!")
            os.mkdir("1000g_vcfs")
            
        # Starting download of required files. 
        if(not os.path.isfile("1000g_vcfs/ALL.chr.bed")
           or not os.path.isfile("1000g_vcfs/ALL.chr.bim")
           or not os.path.isfile("1000g_vcfs/ALL.chr.fam")):
            os.chdir("1000g_vcfs")
            # Downloading LD data. 
            print("LD_SETUP: Initiating download of LD data!")
            os.system(f"aria2c --allow-overwrite=false -x {args.concurrent} -m 0  -i ../urls.txt -d raw")
            
            # Concatenating LD data.
            print("LD_SETUP: Concatenating LD data!")
            os.system(f"bcftools concat -Oz --threads {args.threads} -o 1000g_vcfs/ALL.chr.vcf.gz 1000g_vcfs/raw/ALL.chr*.vcf.gz")
            
            # Removing raw LD data.
            print("LD_SETUP: Removing raw LD data!")
            os.system(f"rm -r raw")
            os.chdir("..")

            # Indexing LD data.
            if (not os.path.isfile("ALL.chr.vcf.gz.tbi")):
                print("LD_SETUP: Indexing LD data!")
                os.system(f"tabix -p vcf 1000g_vcfs/ALL.chr.vcf.gz")
            
            # Reformatting LD data to plink.
            print("LD_SETUP: Converting LD data to plink!")
            os.system(f"plink2 --vcf 1000g_vcfs/ALL.chr.vcf.gz --make-pgen --out 1000g_vcfs/ALL.chr --threads {args.threads}")

            # Removing concatenated data.
            print("LD_SETUP: Removing concatenated LD data!")
            os.system(f"rm -r 1000g_vcfs/ALL.chr.vcf.gz")
            os.system(f"rm -r 1000g_vcfs/ALL.chr.vcf.gz.tbi")

        # Downloading sample information. 
        if (not os.path.isfile("1000g_vcfs/integrated_call_samples_v3.20200731.ALL.ped")):
            os.system(f"wget --timeout=30 --tries=10 --waitretry=5 -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped -O 1000g_vcfs/integrated_call_samples_v3.20200731.ALL.ped")

        # Loading sample information.
        print("LD_SETUP: Loading sample information...")
        sample_info = []
        with open("1000g_vcfs/integrated_call_samples_v3.20200731.ALL.ped", 'r', newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                sample_info.append(row)

        # Annotating sample information.
        print("LD_SETUP: Reformatting sample information...")
        sample_info = pd.DataFrame(sample_info[1:], columns=sample_info[0])

        # Reformat sample information.
        sample_info["FID"] = 0
        sample_info["IID"] = sample_info["Individual ID"]
        sample_info["PAT"] = 0
        sample_info["MAT"] = 0

        # Encode sex information. 
        sample_info["SEX"] = sample_info["Gender"].map({"male": 1, "female": 2}).fillna(0).astype(int)

        # Assign a superpopulation based on population field.
        pop_map = {
            "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
            "AMR": ["MXL", "PUR", "CLM", "PEL"],
            "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
            "EUR": ["CEU", "TSI", "FIN", "GBR", "IBS"],
            "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"]
        }

        # Reverse mapping.
        rev_pop_map = {pop: superpop for superpop, pops in pop_map.items() for pop in pops}
        sample_info["SuperPop"] = sample_info["Population"].map(rev_pop_map)

        # Set phenotype column.
        sample_info["PHENOTYPE"] = -9

        # Select and reorder columns.
        formatted = sample_info[["FID", "IID", "PAT", "MAT", "SEX", "SuperPop", "Population", "PHENOTYPE"]]

        # Delete .ped file.
        os.remove("1000g_vcfs/integrated_call_samples_v3.20200731.ALL.ped")

        # Return formatted sample information.
        return formatted

    def create_pop_lists(formatted: pd.DataFrame):
        """
        Function to create population files.
        Args:
            formatted (pd.DataFrame): Formatted sample information.
        Returns:
            None
        """
        # Validating entered population. 
        pops = ["ALL", "EUR", "AFR", "AMR", "EAS", "SAS",
                      "YRI", "LWK", "GWD", "MSL", "ESN", "ASW",
                      "ACB", "MXL", "PUR", "CLM", "PEL", "CHB",
                      "JPT", "CHS", "CDX", "KHV", "CEU", "TSI",
                      "FIN", "GBR", "IBS", "GIH", "PJL", "BEB",
                      "STU", "ITU"]
        
        # Making a folder for the subsets. 
        if not os.path.isdir("1000g_vcfs/subsets"):
            os.mkdir("1000g_vcfs/subsets")
        os.chdir("1000g_vcfs/subsets")
        subset = None

        # Iterating through populations and creating lists. 
        for pop in pops:
            if (pop in pd.unique(formatted["Population"])):
                subset = formatted[formatted["Population"] == pop]
            elif (pop == "ALL"):
                subset = formatted
            else:
                subset = formatted[formatted["SuperPop"] == pop]
            poplist = subset[["IID"]]
            poplist.to_csv(f"{pop}_list", index = False, header = False, sep = "\t")

def main():
    """
    Main function to handle the execution flow of the pipeline.
    Parses arguments, performs sanity checks, and runs preprocessing and LD data setup.
    """
    # Parse command line arguments.
    arguments = LDcreate.check_args(sys.argv[1:])

    # Create LD data.
    LDcreate.create_pop_lists(LDcreate.create_LD(arguments))

if __name__ == "__main__":
    main()
    print("LD_SETUP: LD data creation complete!")

    