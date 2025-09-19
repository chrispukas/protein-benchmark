import os

from incito_pipeline.util.io import pdb, fasta, yaml_parse, entries

#=====================================
# --- PDB TO PDB PARSE (filtering) ---
#=====================================

class PDBTOPDB():
    def __init__(self, 
                 pdb_file: str = None,
                 chain_map: dict[str, str] = None,
                 chain_filter: list[str] = None):
        self.pdb_file = pdb_file
        self.chain_map = chain_map
        self.chain_filter = chain_filter
        self.records = None
        self.atom_counters = None

        if pdb_file:
            print(f"Unpacking PDB file {pdb_file}")
            r = self._load_pdb(pdb_file=pdb_file)
            print(r)
            self.records, self.atom_counters = r

    def _load_pdb(self,
                 pdb_file: str = None,
                 chain_map: dict[str, str] = None,
                 chain_filter: set[str] = None):
 
        if pdb_file is None:
            pdb_file = self.pdb_file
        if chain_map is None:
            chain_map = self.chain_map
        if chain_filter is None:
            chain_filter = self.chain_filter

        if pdb_file and chain_map:
            print(f"Loading PDB from file {pdb_file}")
            r, _ = pdb.reformat_pdb(*pdb.extract_lines(pdb_file), 
                                            chain_map=chain_map,
                                            chain_filter=chain_filter)
            return r, _

    def save_pdb(self, 
                 output_file: str,
                 chain_map: dict[str, str] = None) -> None:
        if isinstance(chain_map, None):
            chain_map = self.chain_map

        pdb.write_records_to_pdb(records=self.records, 
                                 output_file=output_file)
    def get_sequences_seqres(self) -> dict[str, str]:
        if self.records is None:
            print("No records found!")  
            return None
        return pdb.get_sequences_seqres(self.records)
    
    def get_sequences_records(self) -> dict[str, str]:
        if self.records is None:
            print("No records found!")  
            return None
        return pdb.get_sequences_records(self.records)

    def get_chains(self,
                   records: list = None) -> dict[str, list]:
        if records is None:
            records = self.records

        return pdb.get_chains(records if records is not None else self.records)


#===========================
# --- PDB TO FASTA PARSE ---
#===========================

class PDBTOFASTA():
    def __init__(self, 
                 pdb_file: str,
                 chain_map: dict[str, str] = None,
                 chain_filter: set[str] = None):
        
        self.pdb_file = pdb_file
        self.chain_map = chain_map
        self.chain_filter = chain_filter

    def save_fasta(self,
                pdb_file: str = None,
                output_file: str = None,
                chain_map: dict[str, str] = None,
                chain_filter: set[str] = None):
        if pdb_file is None:
            pdb_file = self.pdb_file
        if output_file is None:
            output_file = self.output_file
        if chain_map is None:
            chain_map = self.chain_map
        if chain_filter is None:
            chain_filter = self.chain_filter

        fasta.fasta_from_pdb(pdb_file,
                             output_file,
                             chain_map=chain_map,
                             chain_filter=chain_filter)
        
#==========================
# --- PDB TO YAML PARSE ---
#==========================

class PDBTOYAML(PDBTOPDB):
    def __init__(self, 
                 pdb_file: str = None, 
                 save_dir: str = None, 
                 chain_map: dict[str, str] = {"A": "A", "H":"H", "L":"L"},
                 msa_path: str = None,
                 templates: list[entries.Template] = [],
                 constraints: list[entries.PocketConstraint] = [],
                 chain_filter: list[str] = None,):
        self.pdb_file = pdb_file
        self.save_dir = save_dir
        self.chain_map = chain_map
        self.msa_path = msa_path
        self.templates = templates
        self.constraints = constraints
        self.chain_filter = chain_filter

        super().__init__(pdb_file=pdb_file,
                         chain_map=chain_map,
                         chain_filter=chain_filter)

    def save_yaml(self,
                  pdb_file: str = None, 
                  save_dir: str = None, 
                  chain_map: dict[str, str] = None,
                  msa_path: str = None,
                  templates: list[entries.Template] = None,
                  constraints: list[entries.PocketConstraint] = None,
                  chain_filter: list[str] = None,):
        if pdb_file is None:
            pdb_file = self.pdb_file
        if save_dir is None:
            save_dir = self.save_dir
        if chain_map is None:
            chain_map = self.chain_map
        if msa_path is None:
            msa_path = self.msa_path
        if templates is None:
            templates = self.templates
        if constraints is None:
            constraints = self.constraints
        if chain_filter is None:
            chain_filter = self.chain_filter

        try:

            os.makedirs(os.path.dirname(save_dir), exist_ok=True)

            yaml_parse.pdb_to_yaml(
                pdb_file=pdb_file,
                save_dir=save_dir,
                chain_map=chain_map,
                msa_path=msa_path,
                templates=templates,
                constraints=constraints,
                chain_filter=chain_filter
            )
        except Exception as e:
            print(f"Error occurred while saving YAML: {e}")

class YAMLPARSE():
    def __init__(self, 
                 records: list = None, 
                 save_dir: str = None, 
                 chain_map: dict[str, str] = {"A": "A", "H":"H", "L":"L"},
                 msa_path: str = None,
                 templates: list[entries.Template] = [],
                 constraints: list[entries.PocketConstraint] = [],
                 chain_filter: list[str] = None,):
        self.records = records
        self.save_dir = save_dir
        self.chain_map = chain_map
        self.msa_path = msa_path
        self.templates = templates
        self.constraints = constraints
        self.chain_filter = chain_filter

    def save_yaml(self,
                  records: list = None, 
                  save_dir: str = None, 
                  chain_map: dict[str, str] = None,
                  msa_path: str = None,
                  templates: list[entries.Template] = None,
                  constraints: list[entries.PocketConstraint] = None,
                  chain_filter: list[str] = None,):
        if records is None:
            records = self.records
        if save_dir is None:
            save_dir = self.save_dir
        if chain_map is None:
            chain_map = self.chain_map
        if msa_path is None:
            msa_path = self.msa_path
        if templates is None:
            templates = self.templates
        if constraints is None:
            constraints = self.constraints
        if chain_filter is None:
            chain_filter = self.chain_filter

        try:

            os.makedirs(os.path.dirname(save_dir), exist_ok=True)

            yaml_parse.save_yaml(
                records=records,
                save_dir=save_dir,
                chain_map=chain_map,
                msa_path=msa_path,
                templates=templates,
                constraints=constraints,
                chain_filter=chain_filter
            )
        except Exception as e:
            print(f"Error occurred while saving YAML: {e}")

    def save_yaml_force_seq(self,
                            sequences: list[str] = None,
                            save_dir: str = None,
                            chain_map: dict[str, str] = None,
                            msa_path: str = None,
                            templates: list[entries.Template] = None,
                            constraints: list[entries.PocketConstraint] = None,
                            chain_filter: list[str] = None,):
        if sequences is None:
            sequences = self.sequences
        if save_dir is None:
            save_dir = self.save_dir
        if chain_map is None:
            chain_map = self.chain_map
        if msa_path is None:
            msa_path = self.msa_path
        if templates is None:
            templates = self.templates
        if constraints is None:
            constraints = self.constraints
        if chain_filter is None:
            chain_filter = self.chain_filter

        try:

            os.makedirs(os.path.dirname(save_dir), exist_ok=True)

            yaml_parse.save_yaml(
                sequences=sequences,
                save_dir=save_dir,
                chain_map=chain_map,
                msa_path=msa_path,
                templates=templates,
                constraints=constraints,
                chain_filter=chain_filter
            )
        except Exception as e:
            print(f"Error occurred while saving YAML: {e}")
