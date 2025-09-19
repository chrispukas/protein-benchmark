
import os
import re
import numpy as np
import copy

from dataclasses import dataclass
from importlib.resources import files

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.vectors import Vector, rotaxis2m, rotmat

from scipy.spatial.transform import Rotation as R

#from incito_pipeline.util import data_utility
#DATASET_PATH = files("incito_pipeline.datasets")

def create_template_pdb(path_to_pdb: str, 
                        output_file: str,
                        chain_ids: list[str]) -> None:
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    pdb_name = os.path.basename(path_to_pdb).split('.')[0]
    pdb = PDBParser(QUIET=True).get_structure(pdb_name, path_to_pdb)

    antigen_id_filter = ["A"]
    antibody_id_filter = ["H", "L"]

    antigen_chains = []
    antibody_chains = []

    # Filter chains based on IDs, using standard antibody/antigen chain IDs
    for chain in pdb.get_chains():
        id = chain.get_id()
        if id in antigen_id_filter:
            antigen_chains.append(chain)
        elif id in antibody_id_filter:
            antibody_chains.append(chain)
    
    # Rotate antigen chains together (for variability)
    if len(antigen_chains) >= 1:
        prand_rotation = np.random.uniform(0, 360, 3)
        prand_translation = np.random.uniform(20, 30, 3)
        rotate_translate_quaternion_proteins(antigen_chains, prand_rotation, prand_translation)

    if len(antibody_chains) >= 1:
        prand_rotation = np.random.uniform(0, 360, 3)
        prand_translation = -np.random.uniform(20, 30, 3)
        rotate_translate_quaternion_proteins(antibody_chains, prand_rotation, prand_translation)

    dump_pdb_from_original(path_to_pdb=path_to_pdb,
                           output_file=output_file,
                           chain_filter=chain_ids)
    return output_file

def dump_pdb_from_original(path_to_pdb: str, 
                           output_file: str, 
                           chain_map: dict[str, str] = {'A':'H', 'B':'L', 'C':'A'},
                           chain_filter: list[int] = None) -> None:
    if not os.path.exists(path_to_pdb) or not os.path.isfile(path_to_pdb):
        print(f"Input PDB file {path_to_pdb} does not exist.")
        return

    original_records = extract_lines(path_to_pdb)
    new_records = reformat_pdb(original_records, 
                               chain_map=chain_map,
                               chain_filter=chain_filter)
    write_records_to_pdb(new_records, output_file)


def reformat_pdb(original_records: list, 
                 chain_map: dict[str, str] = None,
                 chain_filter: list[str] = None,) -> list[any]:
    if chain_filter is not None and len(chain_filter) == 0:
        chain_filter = None

    records_to_copy = ["REMARK", "TITLE", "SEQRES", "TER", "END"]
    reformatted_records = []

    for record in original_records:
        if isinstance(record, AtomLine):
            if chain_filter is None \
                or chain_filter is not None and record.chain_identifier in chain_filter:
                reformatted_records.append(record)

        elif record.record_name in records_to_copy:
            if (isinstance(record, SeqresLine) or isinstance(record, TerminatorLine) or isinstance(record, MissingEntriesLine)) \
                    and (chain_filter is not None and record.chain_identifier not in chain_filter):
                continue

            if isinstance(record, RemarkLine) and record.text.startswith("original_to_current_chain_mapping"):
                record.text = f"original_to_current_chain_mapping={str(chain_map)}"
            reformatted_records.append(record)

    return reformatted_records



def write_records_to_pdb(records: list, 
                         output_file: str):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as out_f:
        for record in records:
            out_f.write(record.original_line + "\n")
    print(f"Wrote {len(records)} records to {output_file}")



# ==============================================================================================
# --- Extract the data from each record in the PDB file, output into readable datastructures ---
# ==============================================================================================

@dataclass
class LineBase():
    original_line: str
    record_name: str

def extract_lines(path_to_pdb: str):
    with open(path_to_pdb) as f:
        lines = f.readlines()
    return format_lines(lines=lines)

def format_lines(lines: list[str]):
    funcs = {"ATOM": get_atom_line_data,
             "SEQRES": get_seqres_line_data,
             "REMARK": get_remark_line_data,
             "TITLE": get_title_line_data,
             "TER": get_terminator_line_data,
             "END": get_default_line_data}
    line_objects = []
    for line in lines:
        for key, func in funcs.items():
            if line.startswith(key):
                data = func(line)
                if data:
                    line_objects.append(data)
                break
    return line_objects

@dataclass
class SeqresLine(LineBase):
    serial_number: int
    chain_identifier: str
    number_of_residues: int
    residue_names: list[str]

def get_seqres_line_data(line: str):
    try:
        extracted_data = {
            "record_name": line[0:6].strip(),
            "serial_number": int(line[8:10].strip()),
            "chain_identifier": line[11].strip(),
            "number_of_residues": int(line[13:17].strip()),
            "residue_names": line[19:].strip().split()
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    record_name = extracted_data["record_name"]
    serial_number = extracted_data["serial_number"]
    chain_identifier = extracted_data["chain_identifier"]
    number_of_residues = extracted_data["number_of_residues"]
    residues = extracted_data["residue_names"]

    return SeqresLine(
        original_line=line.rstrip("\n"),
        record_name=record_name,
        serial_number=serial_number,
        chain_identifier=chain_identifier,
        number_of_residues=number_of_residues,
        residue_names=residues
        )

@dataclass
class TitleLine(LineBase):
    continuation: str
    title_text: str

def get_title_line_data(line: str):
    try:
        extracted_data = {
                "record_name": line[0:6].strip(),
                "continuation": line[8:10].strip(),
                "title_text": line[10:80].strip(),
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    return TitleLine(
        original_line=line.rstrip("\n"),
        record_name=extracted_data["record_name"],
        continuation=extracted_data["continuation"],
        title_text=extracted_data["title_text"]
    )

@dataclass
class AtomLine(LineBase):
    serial_number: int
    atom_name: str
    residue_name: str
    chain_identifier: str
    residue_number: int
    insertion_code: str
    orthagonal_coordinates: tuple[str]
    occupancy: float
    beta_factor: float
    segment_id: str
    element_symbol: str
    charge: str

def get_atom_line_data(line: str):
    try:
        extracted_data = {
            "record_name": line[0:6].strip(),
            "serial_number": int(line[6:11].strip()),
            "atom_name": line[12:16].strip(),
            "residue_name": line[17:20].strip(),
            "chain_identifier": line[21].strip(),
            "residue_number": int(line[22:26].strip()),
            "insertion_code": line[27].strip(),
            "x": float(line[30:38].strip()),
            "y": float(line[38:46].strip()),
            "z": float(line[46:54].strip()),
            "occupancy": float(line[54:60].strip()),
            "beta_factor": float(line[60:66].strip()),
            "segment_id": line[72:76].strip(),
            "element_symbol": line[76:78].strip(),
            "charge": line[78:80].strip(),
        }   
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    record_name: str = extracted_data["record_name"]
    serial_number: int = extracted_data["serial_number"]
    atom_name: str = extracted_data["atom_name"]
    residue_name: str = extracted_data["residue_name"]
    chain_identifier: str = extracted_data["chain_identifier"]
    residue_number: int = extracted_data["residue_number"]
    insertion_code: str = extracted_data["insertion_code"]
    orthagonal_coordinates: tuple[str] = (extracted_data["x"], extracted_data["y"], extracted_data["z"])
    occupancy: float = extracted_data["occupancy"]
    beta_factor: float = extracted_data["beta_factor"]
    segment_id: str = extracted_data["segment_id"]
    element_symbol: str = extracted_data["element_symbol"]
    charge: str = extracted_data["charge"]

    return AtomLine(
        original_line=line.rstrip("\n"),
        record_name=record_name,
        serial_number=serial_number,
        atom_name=atom_name,
        residue_name=residue_name,
        chain_identifier=chain_identifier,
        residue_number=residue_number,
        insertion_code=insertion_code,
        orthagonal_coordinates=orthagonal_coordinates,
        occupancy=occupancy,
        beta_factor=beta_factor,
        segment_id=segment_id,
        element_symbol=element_symbol,
        charge=charge
    )

@dataclass
class RemarkLine(LineBase):
    remark_number: int
    text: str

def get_remark_line_data(line: str):
    try:
        extracted_data = {
            "record_name": line[0:6].strip(),
            "remark_number": int(line[7:10].strip()),
            "text": line[10:].strip()
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None
    
    match extracted_data["remark_number"]: 
        case 465:
            return get_missing_entries_line_data(line)

    return RemarkLine(
        original_line=line.rstrip("\n"),
        record_name=extracted_data["record_name"],
        remark_number=extracted_data["remark_number"],
        text=extracted_data["text"]
    )

@dataclass
class TerminatorLine(LineBase):
    serial_number: int
    residue_name: str
    chain_identifier: str
    residue_number: int
    
def get_terminator_line_data(line: str):
    try:
        extracted_data = {
            "original_line": line.rstrip("\n"),
            "record_name": line[0:6].strip(),
            "serial_number": int(line[6:11].strip()),
            "residue_name": line[17:20].strip(),
            "chain_identifier": line[21].strip(),
            "residue_number": int(line[22:26].strip())
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    return TerminatorLine(
        original_line=extracted_data["original_line"],
        record_name=extracted_data["record_name"],
        serial_number=extracted_data["serial_number"],
        residue_name=extracted_data["residue_name"],
        chain_identifier=extracted_data["chain_identifier"],
        residue_number=extracted_data["residue_number"]
    )

@dataclass
class MissingEntriesLine(LineBase):
    remark_number: int
    residue_name: str
    chain_identifier: str
    res_seq: int

def get_missing_entries_line_data(line: str):
    try:
        extracted_data = {
            "original_line": line.rstrip("\n"),
            "record_name": line[0:6].strip(),
            "remark_number": int(line[7:10].strip()),
            "residue_name": line[18:21].strip(),
            "chain_identifier": line[22].strip(),
            "res_seq": line[23:].strip()
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    return MissingEntriesLine(
        original_line=extracted_data["original_line"],
        record_name=extracted_data["record_name"],
        remark_number=extracted_data["remark_number"],
        residue_name=extracted_data["residue_name"],
        chain_identifier=extracted_data["chain_identifier"],
        res_seq=extracted_data["res_seq"]
    )



@dataclass
class DefaultLine(LineBase):
    text: str

def get_default_line_data(line: str):
    try:
        extracted_data = {
            "original_line": line.rstrip("\n"),
            "record_name": line[0:6].strip(),
            "text": line[7:].strip()
        }
    except ValueError as e:
        print(f"Error parsing DEFAULT line {line}: {e}")
        return None

    return DefaultLine(
        original_line=extracted_data["original_line"],
        record_name=extracted_data["record_name"],
        text=extracted_data["text"]
    )
    

# ====================================================
# --- Geometry helpers for protein transformations ---
# ====================================================


def get_centroid(atoms):
    coords = np.array([atom.coord for atom in atoms])
    return np.mean(coords, axis=0)
def rotate_translate_quaternion_proteins(structures, rotation, translation):
    atoms = [atom for structure in structures for atom in structure.get_atoms()]
    centroid_coordinate = get_centroid(atoms)

    # Method is:
    # 1. Translate atoms to origin
    # 2. Rotate atoms around the origin
    # 3. Translate atoms back to the original centroid position

    # Translate to origin
    transform_atoms(atoms, None, -centroid_coordinate)

    # Rotate atoms around origin
    rotation_matrix = R.from_euler('zyx', rotation, degrees=True).as_matrix()
    transform_atoms(atoms, rotation_matrix, None)

    # Translate back to original centroid position
    transform_atoms(atoms, None, centroid_coordinate + translation)
    return structures
def transform_atoms(atoms, rotation=None, translation=None):
    
    if rotation is None:
        rotation = np.identity(3)
    if translation is None:
        translation = np.zeros(3)
        
    
    for atom in atoms:
        atom.transform(rotation, translation)
    print("Transformed atoms with rotation:\n", rotation, "\nand translation:\n", translation)
