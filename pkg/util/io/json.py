import os
import json

from Bio import SeqIO

def sequence_entry(chain_id: str, seq:str, msa_path: str =None):
    entry = {
        "seq": seq,
        "chain_id": chain_id
    }

    if msa_path is not None:
        entry["msa_path"] = msa_path
    
    return entry


def parse_pdb_to_json(path_to_pdb, json_output_dir):
    # Collect each epitope, and equivalent amino-acid sequences into a list

    pdb_name = os.path.basename(path_to_pdb).split('.')[0]

    records = SeqIO.parse(path_to_pdb, "pdb-seqres")
    entries = []

    for record in records:
        chain_id = str(record.id.split('|')[0])
        seq = str(record.seq)


        entry = sequence_entry(
            chain_id=chain_id,
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

    output_path = os.path.join(json_output_dir, f"{pdb_name}.json")

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)