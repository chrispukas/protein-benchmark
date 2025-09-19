from incito_pipeline.util import data_utility
from importlib.resources import files

from . import data_utility


import os

class BoltzInterface:
    """
    Interface for interacting with BOLTZ, and its output files.
    This class provides methods to convert CIF files to PDB format, manage file paths, and handle querying of BOLTZ.
    """

    def __init__(self, input_cif_dir, output_pdb_dir):
        self.input_cif_dir = input_cif_dir
        self.output_pdb_dir = output_pdb_dir

        # Ensure the output directory exists
        os.makedirs(self.output_pdb_dir, exist_ok=True)
    """
    Boltz interface methods
    """

    def boltz_predict_multiple(self, fasta_file_dir, output_dir, flags=None, test=False):
        """
        Execute a BOLTZ prediction command for multiple FASTA files.
        """
        
        fasta_files = data_utility.get_file_names(fasta_file_dir, filter="combined")
        print(fasta_file_dir)
        if not fasta_files:
            print(f"No FASTA files found in {fasta_file_dir}")
            return
        
        for fasta_file in fasta_files:
            fasta_file_path = os.path.join(fasta_file_dir, fasta_file)
            self.boltz_predict_single(fasta_file_path, output_dir, flags=flags, test=test)


    def boltz_predict_single(self, fasta_file, output_dir, flags=None, test=False):
        """
        Execute a BOLTZ prediction command.
        """
        
        try:
            command = f"boltz predict {fasta_file} --out_dir {output_dir}"
            for flag in flags or []:
                command += f" {flag}"
            print(f"Executing command: {command}")
            
            if not test:
                os.system(command)

        except Exception as e:
            print(f"Error: Failed to execute BOLTZ prediction: {e}")

