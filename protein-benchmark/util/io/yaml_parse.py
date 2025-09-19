import os
import yaml

import incito_pipeline.util.io.entries as e
import incito_pipeline.util.io.pdb as pdb




def sequence_entry(entity_type: str, chain_id: str, sequence: str) -> dict:
    data = {
        "entity_type": entity_type,
        "id": chain_id,
        "sequence": sequence
    }
    return {entity_type: data}


def convert_sequences(sequences: dict[str, str], 
                      msa_path: str) -> list[e.Sequence]:
    seq_entries = []
    for chain_id, sequence in sequences.items():
        sequence = e.Sequence(
            entity_type='protein',
            id=chain_id,
            sequence=sequence,
            msa=msa_path,
            position_modification=None,
        )

        seq_entries.append(sequence)
    return seq_entries


def pdb_to_yaml(pdb_file: str, 
                save_dir: str, 
                chain_map: dict[str, str],
                msa_path: str,
                templates: list[e.Template] = [],
                constraints: list[e.PocketConstraint] = [],
                chain_filter: list[str] = None,
                ) -> None:

    templates = templates or []
    constraints = constraints or []

    records, _ = pdb.reformat_pdb(*pdb.extract_lines(pdb_file),
                                        chain_map=chain_map,
                                        chain_filter=chain_filter)
    save_yaml(records, save_dir, chain_map, msa_path, templates, constraints, chain_filter)

def yaml_dump(data: dict, yaml_file: str, default_flow_style: bool = True) -> None:
    with open(yaml_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=default_flow_style)


def save_yaml(records: list, 
               save_dir: str, 
               chain_map: dict[str, str],
               msa_path: str,
               templates: list[e.Template] = [],
               constraints: list[e.PocketConstraint] = [],
               chain_filter: list[str] = None,
               ) -> None:

    records: list[pdb.LineBase] = records
    sequences: dict[str, str] = pdb.get_sequences_records(records)
    save_yaml_by_sequence(sequences=sequences,
                          save_dir=save_dir,
                          chain_map=chain_map,
                          msa_path=msa_path,
                          templates=templates,
                          constraints=constraints,
                          chain_filter=chain_filter)


def save_yaml_by_sequence(sequences: list[str], 
               save_dir: str, 
               chain_map: dict[str, str],
               msa_path: str,
               templates: list[e.Template] = [],
               constraints: list[e.PocketConstraint] = [],
               chain_filter: list[str] = None,
               ) -> None:

    templates = templates or []
    constraints = constraints or []

    sequences_obj: list[e.Sequence] = convert_sequences(sequences, msa_path)

    sequence_entries = e.get_sequence_entries(sequences_obj)
    template_entries = e.get_template_entries(templates) if templates is not None else []
    constraint_entries = e.get_constraint_entries(constraints) if constraints is not None else []

    data = {
        "sequences": sequence_entries
    }

    if template_entries:
        data["templates"] = template_entries
    if constraint_entries:
        data["constraints"] = constraint_entries

    yaml_dump(data, save_dir)

def yaml_dump(data: dict, yaml_file: str, default_flow_style: bool = True) -> None:
    with open(yaml_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=default_flow_style)