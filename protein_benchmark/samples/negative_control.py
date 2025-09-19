import random
import os
import json

import polars as pf
import numpy as np

import protein_benchmark.extract.pdb_extract as pdb_e
import protein_benchmark.extract.csv_extract as csv_e
import protein_benchmark.io.json as json_p
import protein_benchmark.io.yaml_parse as yaml_p

import protein_benchmark.data_utility as du

from dataclasses import dataclass
from Bio import SeqIO

@dataclass
class Residue:
    sequence_index: int

    residue_number: int

    residue_name_3_letter: str
    residue_name_1_letter: str

def create_negative_control_scramble_json(path_to_pdb_sample, output_dir, cutoff_distance=10):
    # Collect each epitope, and equivalent amino-acid sequences into a list

    pdb_name = os.path.basename(path_to_pdb_sample).split('.')[0]

    residues = pdb_e.get_epitope_residues(
        path_to_pdb=path_to_pdb_sample,
        cutoff_distance=cutoff_distance,
        chain_1_filter=["H", "L"],
        chain_2_filter=["A"]
    )

    records = SeqIO.parse(path_to_pdb_sample, "pdb-seqres")
    print(f"Scrambling {records} records")

    entries = []

    for record in records:
        residue_map = map_residues_from_pdb(record, residues)
        seq = shuffle_sequence_residues(record, residue_map)
        print(seq)

        entry = json_p.sequence_entry(
            chain_id=record.id.split('|')[0],
            seq=seq,
            msa_path=None # No MSA specified
        )

        entries.append(entry)
    
    data = [
        {
            "name": pdb_name,
            "components": entries,
        }
    ]

    output_path = os.path.join(output_dir, f"{pdb_name}_scramble_negative_control.json")

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)






        
def create_negative_control_scramble_yaml(path_to_pdb, output_dir, cutoff_distance=10):
    # Collect each epitope, and equivalent amino-acid sequences into a list

    pdb_name = os.path.basename(path_to_pdb).split('.')[0]

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
        residue_map = map_residues_from_pdb(record, residues)
        seq = shuffle_sequence_residues(record, residue_map)
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

    output_path = os.path.join(output_dir, f"{pdb_name}_scramble_negative_control.yaml")
    yaml_p.yaml_dump(data, output_path)

def shuffle_sequence_residues(record, residue_map):
    sequence_indexes = []
    residue_1_letter = []

    for residue in residue_map:
        sequence_indexes.append(residue.sequence_index)
        residue_1_letter.append(residue.residue_name_1_letter)

    random.shuffle(residue_1_letter)
    seq = list(str(record.seq))

    for i in range(len(sequence_indexes)):
        seq[sequence_indexes[i]] = residue_1_letter[i]
    
    return ''.join(seq)
def map_residues_from_pdb(record, residues) -> list[Residue]:
    print("Attempting to map residues")
    sequences = []

    sequence = record.seq
    chain_id = record.id.split('|')[0]

    for residue in residues:
        if residue.chain_id != chain_id:
            print(f"Residue {residue} does not match chain {chain_id}, skipping.")
            continue
        i = residue.real_residue_number - 1

        residue_name_3_letter = residue.residue_name
        residue_name_1_letter = sequence[i]

        print(f"{residue_name_3_letter} -> {residue_name_1_letter}")

        if not residue_name_3_letter == pdb_e.get_3_letter_code_from_1_letter(residue_name_1_letter):
            print("Error: Residues failed to match")
            continue

        seq = Residue(
            sequence_index = i,
            residue_number = residue.residue_number,
            residue_name_3_letter = residue_name_3_letter,
            residue_name_1_letter = residue_name_1_letter
        )

        sequences.append(seq)

    return sequences


# Negative control samples -> replace chain pairs

def create_negative_control_dataset_replace_single(path_to_pdb: str, 
                                            output_dir: str, 
                                            path_to_id_dataset: str,
                                            path_to_dataset: str, 
                                            chains_to_replace: dict[str: str],
                                            seq_columns: list[str]):

    print("Getting parquet datafile: ", path_to_id_dataset)
    df = csv_e.reformat_db(pf.read_parquet(path_to_id_dataset), seq_columns=seq_columns)
    df_full = pf.read_csv(path_to_dataset)
    
    write_records_to_yaml_nc_replace(path_to_pdb, output_dir, df, df_full, chains_to_replace)


def create_negative_control_dataset_replace_dir(pdb_dirs: list[str], 
                                            output_dir: str, 
                                            path_to_id_dataset: str,
                                            path_to_dataset: str, 
                                            chains_to_replace: dict[str: str],
                                            seq_columns: list[str]) -> None:
    
    print("Getting parquet datafile: ", path_to_id_dataset)
    df = csv_e.reformat_db(pf.read_parquet(path_to_id_dataset), seq_columns=seq_columns)
    df_full = pf.read_csv(path_to_dataset)
    df_full['row_index'] = df_full.index

    for path_to_pdb in pdb_dirs:
        write_records_to_yaml_nc_replace(path_to_pdb, output_dir, df, df_full, chains_to_replace)

    print(f"Recorded {len(pdb_dirs)} records")





def write_records_to_yaml_nc_replace(path_to_pdb: str, 
                                     output_dir: str, 
                                     df, 
                                     df_full, 
                                     chains_to_replace: dict[str: str],):
    pdb_name = os.path.basename(path_to_pdb).split('.')[0]
    print(f"Writing records for negative control chain swap for pdb: {pdb_name}")

    sequence_entries = []
    records = SeqIO.parse(path_to_pdb, "pdb-seqres")
    
    original_antibody_sequences = []

    for record in records:
        chain_id = record.id.split('|')[0]
        seq = str(record.seq)

        if chain_id not in chains_to_replace:
            sequence_entry = yaml_p.sequence_entry(
                entity_type="protein",
                chain_id=chain_id,
                sequence=seq
            )

            sequence_entries.append(sequence_entry)
            continue

        original_antibody_sequences.append((chain_id, seq, len(seq)))
    
    chain_ids, seqs, seq_lengths = zip(*original_antibody_sequences)
    print(f"Passing in lengths: {seq_lengths}")
    target_id = get_id_from_df(df, seq_lengths)

    if target_id is None:
        return

    for chain_id, seq, seq_length in original_antibody_sequences:
        new_seq = get_seq_from_id(df_full, target_id, chains_to_replace[chain_id])

        print(f"Selected target id: {target_id}")
        print(f"Replaced {chain_id} sequences: {str(seq)} -> {new_seq}")
        print(f"Original chain length: {len(seq)}, new chain length: {len(new_seq)}, difference: {len(seq) - len(new_seq)}\n")

        sequence_entry = yaml_p.sequence_entry(
            entity_type="protein",
            chain_id=chain_id,
            sequence=new_seq
        )

        sequence_entries.append(sequence_entry)

    data = {
        "sequences": sequence_entries,
    }

    output_path = os.path.join(output_dir, f"{pdb_name}_replace_negative_control.yaml")
    yaml_p.yaml_dump(data, output_path)


def get_id_from_df(df, sequence_lengths: list[int]):
    try:
        query = df.loc[tuple(sequence_lengths), 'row_index']
        print(type(query))
        print(query)
        if isinstance(query, list):
            ids = query
        else:
            ids = query.to_pandas().to_list()[0]
        if (ids is None) or (len(ids) == 0):
            print(f"No ids found for lengths: {sequence_lengths}")
            return None
        rand_id = np.random.choice(ids)
        return rand_id
    except KeyError:
        print(f"KeyError: {tuple(sequence_lengths)} not found in DataFrame")
        return None


def get_seq_from_id(df, id, type_input):
    query = df.loc[df['row_index'] == id, type_input]
    print(type(query))
    if isinstance(query, list):
        val = query[0]
    else:
        val = query.to_pandas().to_list()[0]
    return str(val)
