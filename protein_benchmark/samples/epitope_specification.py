import os
import numpy as np
import random

import protein_benchmark.extract.pdb_extract as pdb_extract
import protein_benchmark.samples.alanine_swap as alanine_swap
import protein_benchmark.io.io as io
import protein_benchmark.io.entries as e



def format_epitope_residues(residues: tuple[pdb_extract.epitope_residue]) -> set:
    formatted = set()
    for res in residues:
        entry = [
            str(res.chain_id),
            int(res.real_residue_number) # Position in fasta file
        ]
        formatted.add(tuple(entry))
    return formatted



def create_pockets(pdb_file: str,
                   epitope_cutoff: float = 10.0,
                   pocket_constraint: float = 5.0,
                   enforce_pocket: bool = False,
                   chain_filter_model: list[str] = None,
                   chain_filter_native: list[str] = None,
                   chain_map: dict[str, str] = None,
                   percentage: float = 1.0
                   ) -> list[e.PocketConstraint]:
    
    pdb2pdb = io.PDBTOPDB(pdb_file=pdb_file,
                          chain_filter=["A", "H", "L"],
                          chain_map=chain_map)
    records = pdb2pdb.records

    epitope_residues = alanine_swap.get_contact_points(pdb2pdb.get_chains(records=records),
                                                       chain_1_filter=chain_filter_model,
                                                       chain_2_filter=chain_filter_native,
                                                       cutoff_distance=epitope_cutoff)
    mapping: dict[str, set[tuple[str, int]]] = {}
    for chain_id_1, fasta_index_1, _, chain_id_2, _, _ in epitope_residues:
        mapping.setdefault(chain_id_2, set()).add((chain_id_1, fasta_index_1))

    # Flatten all residues across chains
    all_contacts = [(chain_id, contact) for chain_id, contacts in mapping.items() for contact in contacts]

    # Number to keep overall
    n_keep = max(1, int(len(all_contacts) * percentage))

    # Randomly sample from the global pool
    sampled = random.sample(all_contacts, k=n_keep)

    # Rebuild mapping only with sampled residues
    sampled_mapping: dict[str, set[tuple[str, int]]] = {}
    for chain_id, contact in sampled:
        sampled_mapping.setdefault(chain_id, set()).add(contact)

    pockets: list[e.PocketConstraint] = []
    for chain_id, contacts in sampled_mapping.items():
        print(
            f"Chain {chain_id}: kept {len(contacts)}/{len(mapping[chain_id])} residues "
            f"({(len(contacts)/len(mapping[chain_id]))*100:.2f}%)"
        )

        pocket = e.PocketConstraint(
            binder=str(chain_id),
            contacts=contacts,
            force=enforce_pocket,
            max_distance=pocket_constraint
        )
        pockets.append(pocket)

    return pockets


def create_templates(path_to_templates: list[str],
                     force: bool = None,
                     threshold: float = None,
                     chain_ids: list[list[str]] = None,
                     template_ids: list[list[str]] = None,):
    templates: list[e.Template] = []
    path_to_templates = path_to_templates or []
    for i, template_path in enumerate(path_to_templates):

        if template_path.endswith('.pdb'):
            template = e.Template(
                pdb_path=template_path,
                force=force,
                threshold=threshold,
                chain_id=chain_ids[i],
                template_id=template_ids[i]
            )
        elif template_path.endswith('.cif'):
            template = e.Template(
                cif_path=template_path,
                force=force,
                threshold=threshold,
                chain_id=chain_ids[i],
                template_id=template_ids[i]
            )
        templates.append(template)
    return templates


def create_yaml_with_retrieved_epitopes(pdb_file: str, 
                                        save_dir: str, 
                                        epitope_cutoff: float = 10.0,
                                        pocket_constraint_threshold: float = None,
                                        template_constraint_threshold: float = None,
                                        enforce_pocket: bool = None,
                                        msa: str = None,
                                        template_paths: list[str] = None, 
                                        enforce_templates: bool = None,
                                        chain_ids: list[list[str]] = None,
                                        template_ids: list[list[str]] = None,
                                        chain_filter_model: list[str] = None,
                                        chain_filter_native: list[str] = None,
                                        chain_map: dict[str, str] = None,
                                        pocket_percentage: float = 1.0,
                                        ) -> None:
    os.makedirs(os.path.dirname(save_dir), exist_ok=True)
    p2y = io.PDBTOYAML()

    if enforce_pocket and pocket_constraint_threshold:
        pockets = create_pockets(pdb_file,
                                epitope_cutoff=epitope_cutoff,
                                pocket_constraint=pocket_constraint_threshold,
                                enforce_pocket=enforce_pocket,
                                chain_filter_model=chain_filter_model,
                                chain_filter_native=chain_filter_native,
                                chain_map=chain_map,
                                percentage=pocket_percentage)
    else:
        pockets = None

    templates = create_templates(path_to_templates=template_paths,
                                  force=enforce_templates,
                                  threshold=template_constraint_threshold,
                                  chain_ids=chain_ids,
                                  template_ids=template_ids)
    
    output_path = os.path.join(save_dir, f"{os.path.basename(pdb_file).split('.')[0]}.yaml")
    p2y.save_yaml(pdb_file=pdb_file,
                  save_dir=output_path,
                  constraints=pockets,
                  templates=templates,
                  msa_path=msa,)



def create_yaml_with_retrieved_epitopes_multiple(pdb_files: list[str], 
                                                 save_dir: str,
                                                 cutoff_distance: float = 10.0,
                                                 enforce_pocket: bool = False,
                                                 enforce_template: bool = False) -> None:
    for pdb_file in pdb_files:
        create_yaml_with_retrieved_epitopes(pdb_file, 
                                            save_dir, 
                                            cutoff_distance=cutoff_distance, 
                                            enforce_pocket=enforce_pocket,
                                            enforce_template=enforce_template)
                                            
    print(f"{len(pdb_files)} YAML files created in {save_dir}")