import configparser

class Constants:
    """
    A class of constants (mainly to store important directories/inputs).
    """
    def __init__(self, test= "False"): 
        """
        Instantiation function which reads all of the config options into python.
        """
        self._config = configparser.ConfigParser()
        if test == "True" or test == "true" or test == "T" or test == "t":
            self._config.read("config_test.ini")
        else:    
            self._config.read("config.ini")

        # Each stores the corresponding config information.
        self._GENES = self._config.get("genes", "genes")
        self._SEQID_DIR = self._config.get("seqIDdir", "seqIDdir")
        self._DATA_DIR = self._config.get("data_dir", "data_dir")
        self._GWAS_DIR = self._config.get("GWASdir", "GWASdir")
        self._CHR_INPUT = self._config.get("CHR_input", "CHR_input")
        self._BP_INPUT = self._config.get("BP_input", "BP_input")
        self._A1_INPUT = self._config.get("A1_input", "A1_input")
        self._A2_INPUT = self._config.get("A2_input", "A2_input")
        self._BETA_INPUT = self._config.get("BETA_input", "BETA_input")
        self._SE_INPUT = self._config.get("SE_input", "SE_input")
        
    # A bunch of getters for the different constants. 
    def genes(self): return self._GENES
    def seqID_dir(self): return self._SEQID_DIR
    def data_dir(self): return self._DATA_DIR
    def GWASdir(self): return self._GWAS_DIR
    def CHR_input(self): return self._CHR_INPUT
    def BP_input(self): return self._BP_INPUT
    def A1_input(self): return self._A1_INPUT
    def A2_input(self): return self._A2_INPUT
    def BETA_input(self): return self._BETA_INPUT
    def SE_input(self): return self._SE_INPUT