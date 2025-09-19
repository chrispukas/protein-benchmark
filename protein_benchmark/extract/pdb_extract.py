import os

from dataclasses import dataclass

from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO


@dataclass(frozen=True)
class c_alpha:
    residue_number: int
    real_residue_number: int
    residue_name: str
    
    coordinate: tuple[float, float, float]

@dataclass(frozen=True)
class chain_c_alpha:
    chain_id: str
    c_alpha_atoms: tuple[c_alpha, ...]

@dataclass(frozen=True)
class epitope_residue:    
    chain_id: str

    residue_number: int
    real_residue_number: int

    residue_name: str

def get_epitope_residues(
        path_to_pdb: str, 
        cutoff_distance: float, 
        chain_1_filter: list[str] = ["A"], 
        chain_2_filter: list[str] = ["H", "L"],
        ) -> tuple[epitope_residue]:

    pdb = get_pdb_structure(path_to_pdb)
    chain_1, chain_2 = sort_chains(pdb, chain_1_filter=chain_1_filter, chain_2_filter=chain_2_filter)

    epitope_residues = []
    logged_residues = []

    # Iterate through each antigen chain
    for chain_1_item in chain_1:
        chain_1_id = chain_1_item.chain_id
        for chain_1_atom in chain_1_item.c_alpha_atoms:

            # Iterate through each antigen chain
            for chain_2_item in chain_2:
                for chain_2_atom in chain_2_item.c_alpha_atoms:
                    residue_number = chain_2_atom.residue_number
                    distance_between_ca = get_sq_distance_between_atoms(chain_1_atom, chain_2_atom)
                    if distance_between_ca < cutoff_distance and residue_number not in logged_residues:
                        residue = epitope_residue(
                            chain_id=chain_1_id,
                            residue_number=residue_number,
                            real_residue_number=chain_1_atom.real_residue_number,
                            residue_name=chain_1_atom.residue_name,
                        )
                        epitope_residues.append(residue)
                        logged_residues.append(residue_number)

    print(f"Found {len(epitope_residues)} epitope residues in the antibody/antigen chains.")
    return tuple(epitope_residues)


def get_sq_distance_between_atoms(atom1: c_alpha, atom2: c_alpha) -> float:
    coord_1 = atom1.coordinate
    coord_2 = atom2.coordinate

    a_sq = (coord_1[0] - coord_2[0]) ** 2
    b_sq = (coord_1[1] - coord_2[1]) ** 2
    c_sq = (coord_1[2] - coord_2[2]) ** 2

    return (a_sq + b_sq + c_sq) ** 0.5



def sort_chains(pdb_structure, 
                chain_1_filter: list[str], 
                chain_2_filter: list[str],
                skip_empty_residues: bool = False):

    print(f"Sorting chains with filters: {chain_1_filter} and {chain_2_filter}")

    chain_1 = [] # chain id: H, L
    chain_2 = [] # chain id: A

    for chain in pdb_structure.get_chains():
        chain_id = chain.get_id()

        c_alpha_atoms = extract_c_alpha_coordinates_from_chain(chain)

        c_alpha_chain = chain_c_alpha(
            chain_id=chain_id,
            c_alpha_atoms=c_alpha_atoms,
        )

        if chain_id in chain_1_filter:
            chain_1.append(c_alpha_chain)
        elif chain_id in chain_2_filter:
            chain_2.append(c_alpha_chain)
        else:
            print(f"Chain {chain_id} not in filter, skipping.")
    return chain_1, chain_2

def extract_c_alpha_coordinates_from_chain(chain) -> tuple[c_alpha]:
    c_alpha_atoms = []

    residues = list(chain.get_residues())
    offset = residue_ids(residues)[0]
    current_residue = 1

    for residue in residues:
            
            atoms = list(residue.get_atoms())
            residue_id = residue.get_id()
            residue_number = residue_id[1]
            for atom in atoms:
                id = atom.get_id()
                if id.upper() == "CA":

                    coord = atom.get_coord()
                    c_a = c_alpha(
                        residue_name=residue.get_resname(),
                        residue_number=residue_number,
                        real_residue_number=residue_number - offset + 1,
                        coordinate=coord,
                    )
                    c_alpha_atoms.append(c_a)
            current_residue += 1
    return tuple(c_alpha_atoms)

def residue_ids(residues):
    return sorted(res.get_id()[1] for res in residues)

def get_pdb_structure(path_to_pdb: str):
    pdb_name = os.path.basename(path_to_pdb)
    structure = PDBParser().get_structure(pdb_name, path_to_pdb)
    return structure



@dataclass(frozen=True)
class sequence:
    sequence_entry: str
    chain_id: str
    sequence: str


def get_sequence_entries(path_to_pdb: str, 
                         chain_ids: tuple[str] = None) -> list:
    sequences = []

    parser = SeqIO.parse(path_to_pdb, "pdb-seqres")

    for record in parser:
        chain_id = record.id.split('|')[0]
        seq = record.seq

        if chain_ids is not None and chain_id not in chain_ids:
            continue

        seq_obj = sequence(
            sequence_entry="protein",
            chain_id=chain_id,
            sequence=seq,
        )

        # Add sequence entry
        sequences.append(seq_obj)

    return sequences


def get_3_letter_code_from_1_letter(code: str) -> str:
    mapping = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
        "C": "CYS", "E": "GLU", "Q": "GLN", "G": "GLY",
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    return mapping.get(code.upper(), "UNK")  # "UNK" for unknown