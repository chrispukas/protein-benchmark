from incito_pipeline.util.data_utility import get_dirs, get_file_names

import os
import copy
import gemmi
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Structure, Model, MMCIFIO, MMCIF2Dict


def convert_pdb_from_cif_nested(boltz_output_dir, output_dir):
    prediction_dirs = get_dirs(boltz_output_dir, filter="_combined")

    for prediction_dir in prediction_dirs:
        pdb_file_name = prediction_dir.split('_combined')[0]
        prediction_path_cif = os.path.join(boltz_output_dir, prediction_dir)
        output_pdb_dir = os.path.join(output_dir, f"{pdb_file_name}")
        
        if not os.path.exists(output_pdb_dir):
            os.makedirs(output_pdb_dir)

        cif_files = get_file_names(prediction_path_cif, ".cif")
        for cif_file in cif_files:
            path_to_cif = os.path.join(prediction_path_cif, cif_file)
            path_to_pdb_output = os.path.join(output_pdb_dir, f"{cif_file.split('.')[0]}.pdb")
            convert_pdb_from_cif(path_to_cif, path_to_pdb_output)
            print("Converted CIF to PDB for:", cif_file)

        print(f"Converted {prediction_dir} to PDB format in {output_pdb_dir}")


def convert_pdb_from_cif(cif_file, output_file_path, override: bool = False, simplify_name: bool = True):
    basename = os.path.basename(output_file_path)
    dirname = os.path.dirname(output_file_path)
    sep = basename.split('_')
    name = f"{sep[0]}_{sep[1]}_{sep[2].split('.')[0]}.pdb"

    output_file_path = os.path.join(dirname, name) if simplify_name else basename
    
    if os.path.exists(output_file_path) and not override:
        print("Skipping, file already exists.")
        return

    pdb_parse(cif_file, output_file_path)

def pdb_parse(cif_file, pdb_file, override: bool = False, simplify_name: bool = True):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(cif_file).split('.')[0], cif_file)

    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_file)

    print(f"Converted {cif_file} to {pdb_file}")


def convert_pdb_from_cif_dataset(input_cif_dir: str, output_pdb_dir: str, override: bool = False, simplify_name: bool = True):
    # Get all directories in the BOLTZ output directory
    boltz_output_cif_names = get_dirs(input_cif_dir)
    sucess_count = 0

    for dir in boltz_output_cif_names:
        predictions_dir = os.path.join(input_cif_dir, dir)
        cif_file_names = get_file_names(predictions_dir, filter=".cif")

        for cif_file in cif_file_names:
            cif_file_path = os.path.join(predictions_dir, cif_file)
            new_file_name = f"{dir}"

            save_dir = os.path.join(output_pdb_dir, new_file_name)

            convert_pdb_from_cif(cif_file_path, f"{save_dir.split('_combined')[0]}.pdb", override=override, simplify_name=simplify_name)
            print(f"  Converted {cif_file} to PDB format and saved to {save_dir}")
            sucess_count += 1
    print(f"# Converted {sucess_count} cif files to pdb successfully.")





import incito_pipeline.util.samples.templates as templates

def separate_pdb(path_to_pdb: str, 
                 output_dir: str):
    ab_chain_ids = ["H", "L"]
    ag_chain_ids = ["A"]

    templates.dump_pdb_from_original(path_to_pdb=path_to_pdb,
                                      output_file=os.path.join(output_dir, f"ab.pdb"),
                                      chain_map={'H':'A', 'L':'B'},
                                      chain_filter=ab_chain_ids)
    templates.dump_pdb_from_original(path_to_pdb=path_to_pdb,
                                      output_file=os.path.join(output_dir, f"ag.pdb"),
                                      chain_map={'A':'C'},
                                      chain_filter=ag_chain_ids)






from dataclasses import dataclass



def reformat_pdb(original_records: list,
                 atom_counters,
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

    return reformatted_records, atom_counters



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
    atom_counters = {}
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
                if key == "ATOM":
                    
                    data, atom_counters = safe_parse(func, line=line, atom_counters=atom_counters)
                else:
                    data = safe_parse(func, line=line)

                if data:
                    line_objects.append(data)
                break
    return line_objects, atom_counters



def safe_parse(func, **kwargs):
    try:
        return func(**kwargs)
    except Exception as e:
        print(f"Error parsing line with {func.__name__}: {e}")
        return None


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
    fasta_index: int
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

def get_atom_line_data(line: str, atom_counters):
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

    chain_id_upper = chain_identifier.upper()
    key = f"{residue_number}{insertion_code}"

    if chain_id_upper not in atom_counters:
        atom_counters[chain_id_upper] = {}

    if key not in atom_counters[chain_id_upper]:
        atom_counters[chain_id_upper][key] = len(atom_counters[chain_id_upper])

    fasta_index = atom_counters[chain_id_upper][key]


    return AtomLine(
        fasta_index=fasta_index,
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
    ), atom_counters

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


# ============================
#   --- Helper Functions ---
# ============================




def get_chains(records: list) -> dict[str, list]:
    hmap = {}
    for record in records:
        if not isinstance(record, AtomLine):
            continue
        
        chain_id = record.chain_identifier
        if chain_id not in hmap:
            hmap[chain_id] = []

        hmap[chain_id].append(record)

    return hmap


def get_sequences_seqres(records: list) -> dict[str, str]:
    sequences = {}
    for record in records:
        if not isinstance(record, SeqresLine):
            continue
        chain_id = record.chain_identifier
        residues = record.residue_names
        for res in residues:
            res_1 = _3_to_1(res)
            if chain_id not in sequences:
                sequences[chain_id] = []
            sequences[chain_id].append(res_1)

    for chain_id, residues in sequences.items():
        sequences[chain_id] = "".join(residues)
        
    return sequences


def get_sequences_records(records: list) -> dict[str, str]:
    seen: dict[str, set[int]] = {}
    sequences = {}

    for record in records:
        if not isinstance(record, AtomLine):
            continue
        chain_id = record.chain_identifier
        r_num = record.residue_number
        r_name = record.residue_name
        i_code = record.insertion_code

        key = f"{r_num}{i_code}"

        if chain_id not in seen:
            seen[chain_id] = set()
        if key not in seen[chain_id]:
            if chain_id not in sequences:
                sequences[chain_id] = ""
            sequences[chain_id] += _3_to_1(r_name)
        seen[chain_id].add(key)
    return sequences


def _1_to_3(code: str) -> str:
    mapping = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
        "C": "CYS", "E": "GLU", "Q": "GLN", "G": "GLY",
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    return mapping.get(code.upper(), "UNK")

def _3_to_1(code3: str) -> str:
    mapping = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
        "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    return mapping.get(code3.upper(), "X")




