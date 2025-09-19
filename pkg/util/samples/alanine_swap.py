import incito_pipeline.util.io.io as io
import incito_pipeline.util.io.yaml_parse as yaml_p
import incito_pipeline.util.io.entries as e
import incito_pipeline.util.io.pdb as pdb_p
import incito_pipeline.util.extract.pdb_extract as pdb_e
import incito_pipeline.util.samples.negative_control as nc

import copy
import scipy
import os

from Bio import SeqIO



def create_negative_control_scramble_yaml(path_to_pdb, output_dir, cutoff_distance=10):
    # Collect each epitope, and equivalent amino-acid sequences into a list
    residues = pdb_e.get_epitope_residues(
        path_to_pdb=path_to_pdb,
        cutoff_distance=cutoff_distance,
        chain_1_filter=["H", "L"],
        chain_2_filter=["A"]
    )

    records = SeqIO.parse(path_to_pdb, "pdb-seqres")
    print(f"Scrambling {records} records")

    sequence_entries = []

    for record in records:
        residue_map = nc.map_residues_from_pdb(record, residues)
        seq = replace_sequence_residues(record, residue_map)
        chain_id = record.id.split('|')[0]

        sequence_entry = yaml_p.sequence_entry(
            entity_type="protein",
            chain_id=chain_id,
            sequence=seq
        )

        sequence_entries.append(sequence_entry)
    
    data = {
        "sequences": sequence_entries,
    }
    yaml_p.yaml_dump(data, output_dir)


def replace_sequence_residues(record, residue_map):
    sequence_indexes = []

    for residue in residue_map:
        sequence_indexes.append(residue.sequence_index)

    seq = list(str(record.seq))

    for i in range(len(sequence_indexes)):
        seq[sequence_indexes[i]] = "A"
    
    return ''.join(seq)






def replace_residues_on_epitope(pdb_path: str,
                     save_file: str,
                     replace_with: str = "A",
                     chain_map: dict[str, str] = None,
                     epitope_cutoff_distance: float = 10,
                     chain_1_filter: list[str] = ["A"],
                     chain_2_filter: list[str] = ["H", "L"],) -> None:
    pdb2pdb: io.PDBTOPDB = io.PDBTOPDB(pdb_path, chain_map=chain_map)
    records: list = pdb2pdb.records
    sequences: dict[str, str] = pdb2pdb.get_sequences_records()

    chains: dict[str, list] = pdb2pdb.get_chains(records=records)
    contact_points = get_contact_points(chains, 
                                        chain_1_filter=chain_1_filter,
                                        chain_2_filter=chain_2_filter,
                                        cutoff_distance=epitope_cutoff_distance)
    print(f"   {sequences["A"]}")
    l = list(sequences["A"])
    seen_ids = set()
    for chain_id, residue_number_1, residue_id, _, _, _ in contact_points:
        if chain_id == "A" and residue_number_1 not in seen_ids:
            orig = l[residue_number_1]
            if orig == residue_id:
                l[residue_number_1] = replace_with
            else:
                print(f" Error: Residue mismatch {orig} != {residue_id}, quitting!")
            seen_ids.add(residue_number_1)
    sequences["A"] = ''.join(l)
    print(f"   {sequences["A"]}")

    pdb2yaml: io.YAMLPARSE = io.YAMLPARSE()
    pdb2yaml.save_yaml(
        records=records,
        save_dir=save_file,
        chain_map=chain_map,
        msa_path=None,
        templates=None,
        constraints=None,
        chain_filter=None
    )


import random

def scramble_residues_on_epitope(pdb_path: str,
                                 save_file: str,
                                 chain_map: dict[str, str] = None,
                                 epitope_cutoff_distance: float = 10,
                                 chain_1_filter: list[str] = ["A"],
                                 chain_2_filter: list[str] = ["H", "L"]) -> None:
    """
    Scramble (shuffle) residues at the epitope of chain 'A' that are
    within epitope_cutoff_distance to chains in chain_2_filter.
    Saves the modified PDB/YAML to save_file.
    """
    
    # Load PDB and extract sequences & records
    pdb2pdb: io.PDBTOPDB = io.PDBTOPDB(pdb_path, chain_map=chain_map)
    records: list = pdb2pdb.records
    sequences: dict[str, str] = pdb2pdb.get_sequences_records()

    # Extract chains
    chains: dict[str, list] = pdb2pdb.get_chains(records=records)

    # Get epitope contact points
    contact_points = get_contact_points(
        chains, 
        chain_1_filter=chain_1_filter,
        chain_2_filter=chain_2_filter,
        cutoff_distance=epitope_cutoff_distance
    )

    print(f"Original sequence (chain A): {sequences['A']}")

    # Collect epitope indices and residues
    contact_indices = []
    contact_residues = []
    seen_ids = set()
    for chain_id, residue_number_1, residue_id, _, _, _ in contact_points:
        if chain_id == "A" and residue_number_1 not in seen_ids:
            if sequences["A"][residue_number_1] != residue_id:
                print(f"Error: Residue mismatch {sequences['A'][residue_number_1]} != {residue_id}, quitting!")
                return
            contact_indices.append(residue_number_1)
            contact_residues.append(residue_id)
            seen_ids.add(residue_number_1)

    # Shuffle the residues
    random.shuffle(contact_residues)

    # Update the sequence
    sequence_list = list(sequences["A"])
    for idx, new_residue in zip(contact_indices, contact_residues):
        sequence_list[idx] = new_residue
    sequences["A"] = ''.join(sequence_list)
    print(f"Shuffled sequence (chain A): {sequences['A']}")

    # Update the PDB records to match shuffled residues
    for record in records:
        if record.get("chain_id") == "A" and record.get("residue_number") in contact_indices:
            idx = contact_indices.index(record["residue_number"])
            record["resName"] = contact_residues[idx]  # assuming single-letter to three-letter mapping handled elsewhere

    # Save as YAML
    pdb2yaml: io.YAMLPARSE = io.YAMLPARSE()
    pdb2yaml.save_yaml(
        records=records,
        save_dir=save_file,
        chain_map=chain_map,
        msa_path=None,
        templates=None,
        constraints=None,
        chain_filter=None
    )





# Implement KD-tree for n log n complexity
#  Set datatype to be float32
#  K-d tree with O(nlog n) complexity

def get_contact_points(chains: dict[str, list],
                       chain_1_filter: list[str] = ["A"],
                       chain_2_filter: list[str] = ["H", "L"],
                       cutoff_distance: float = 10) -> list[tuple[str, int, str, str, int, str]]:
    contact_points = []

    for i, (chain_id_1, chain_1_atoms) in enumerate(chains.items()):
        if chain_id_1 not in chain_1_filter:
            continue
        for j, (chain_id_2, chain_2_atoms) in enumerate(chains.items()):
            if chain_id_2 not in chain_2_filter:
                continue

            for atom_1 in chain_1_atoms:
                if atom_1.atom_name == "":
                    continue
                for atom_2 in chain_2_atoms:
                    if atom_2.atom_name == "":
                        continue
                    d = distance(atom_1, atom_2)
                    if d <= cutoff_distance:
                        contact_points.append((
                            chain_id_1, atom_1.fasta_index, pdb_p._3_to_1(atom_1.residue_name),
                            chain_id_2, atom_2.fasta_index, pdb_p._3_to_1(atom_2.residue_name)
                        ))

    return list(contact_points)

def get_contact_points_kd(chains: dict[str, list],
                        chain_1_filter: list[str] = ["A"],
                        chain_2_filter: list[str] = ["H", "L"],
                        cutoff_distance: float = 10) -> list:
    chains[chain_1_filter[0]]
    

def distance(atom_1: pdb_p.AtomLine, 
             atom_2: pdb_p.AtomLine) -> float:
    
    c_1 = atom_1.orthagonal_coordinates
    c_2 = atom_2.orthagonal_coordinates

    return ((c_1[0] - c_2[0])**2 + (c_1[1] - c_2[1])**2 + (c_1[2] - c_2[2])**2)**0.5 # Euler coordinates.