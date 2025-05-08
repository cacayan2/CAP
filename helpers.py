from imports import *

class Helpers:
    """
    A class of Helper functions.
    """
    @abstractmethod
    def check_arg(args = None):
        """
        Function to parse command line arguments (implemented).
        Args:
            args (list, optional): List of arguments passed from the command line.
        Returns: 
            argpase.namespace: Parsed arguments object.
        """
        parser = argparse.ArgumentParser(description="wrapper.py")
        parser.add_argument("-p", "--process", default = "pqtl",
                            help = "If using eqtl data, add process flag, otherwise pqtl is the default.",
                            required = False)
        parser.add_argument("-t", "--test",
                            help = "Test the pipeline with test data.",
                            required = False)
        parser.add_argument("-o", "--output",
                            help = "Name of directory where you want output, this will be made for you.",
                            required = True)
        parser.add_argument("-p", "--population",
                            help = "Name of population to pull LD data from", required = True)
        parser.add_argument("-l", "--lddir",
                            help = "Location of LD data.", required = True)
        parser.add_argument("-c", "--threads",
                            help = "Number of threads to use.", required = "True")
        return parser.parse_args(args)
    
    @abstractmethod
    def sanity_check(arguments: argparse.Namespace):
        """
        Function to perform sanity checks on command line arguments (implemented).
        Args: 
            arguments (argparse.Namespace): Parsed arguments object.
        Returns:
            None
        """
        
        # List of valid command line arguments.
        valid_process = ["eqtl", "pqtl"]
        valid_test = ["True", "true", "False", "false", "T", "t", "F", "f"]
        valid_population = ["ALL", "EUR", "AFR", "AMR", "EAS", "SAS", 
                            "YRI", "LWK", "GWD", "MSL", "ESN", "ASW",
                            "ACB", "MXL", "PUR", "CLM", "PEL", "GIH", 
                            "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", 
                            "TSI", "FIN", "IBS", "GBR", "IBS", "GIH",
                            "PJL", "BEB", "STU", "ITU"]

        # Validating arguments.
        if arguments.process not in valid_process:
            raise SystemExit("WRAPPER: Process flag must be either 'eqtl' or 'pqtl'.")
        if arguments.test not in valid_test:
            raise SystemExit("WRAPPER: Test flag must be True or False.")
        if arguments.population not in valid_population:
            raise SystemExit("WRAPPER: Population flag must be one of those listed in the README.")
        if not arguments.threads.isdigit():
            raise SystemExit("WRAPPER: Threads flag must be an integer.")

        # Printing to output.
        print(f"WRAPPER: Running process wrapper.py -p {arguments.process} -t {arguments.test} -o {arguments.output} -s {arguments.population} -c {arguments.threads}")
    
    @abstractmethod
    def refresh_output_directory(output_dir: str):
        """
        Function to remove existing output directory and create a new one.
        """
        if os.path.exists(output_dir):
            os.system(f"rm -r {output_dir}")
        os.system(f"mkdir {output_dir}")
