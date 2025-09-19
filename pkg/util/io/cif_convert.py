from __future__ import annotations

import argparse
import multiprocessing
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from datetime import date
from functools import partial
from pathlib import Path
from typing import Iterable, List

import gemmi
import yaml
from Bio import Align, SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import (MMCIFIO, MMCIFParser, NeighborSearch, PDBParser)
from Bio.PDB.PDBIO import StructureIO
from Bio.PDB.Model import Model as BioModel
from Bio.PDB.Structure import Structure as BioStructure
from Bio.SeqUtils import seq1

# --- Custom Data Classes for YAML output ---
class ContactList(list): pass
class InlineList(list): pass
class ContactPair(list): pass

# --- Scenario Definition ---
@dataclass
class Scenario:
    """Defines a single processing scenario with boolean flags."""
    name: str
    use_msa: bool
    use_template: bool
    use_pocket: bool
    output_dir: Path = field(init=False)


# ==============================================================================
# --- Helper Functions (Global Scope) ---
# ==============================================================================

def is_valid_chain(chain_id) -> bool:
    """Defensively checks if a chain ID is valid, handling various null types."""
    if chain_id is None:
        return False
    
    if hasattr(chain_id, 'item'):
        chain_id = chain_id.item() if getattr(chain_id, 'size', 1) == 1 else None

    if chain_id is None:
        return False
    
    chain_str = str(chain_id).strip()
    return bool(chain_str) and chain_str.upper() not in ['NA', 'NAN', '']

def load_chain_mapping(mapping_file_path: Path) -> dict:
    """Loads and processes the chain mapping CSV file using Polars."""
    import polars as pl
    try:
        df = pl.read_csv(mapping_file_path, null_values=[''])
        rename_dict = {'pdb ID': 'pdb_id', 'A chain(s)': 'antigen', 'H chain(s)': 'heavy', 'L chain(s)': 'light'}
        existing_cols = [col for col in rename_dict if col in df.columns]
        df = df.select(existing_cols).rename({k: v for k, v in rename_dict.items() if k in existing_cols})

        expressions = [pl.col('pdb_id').str.to_lowercase()]
        for col_name in ['antigen', 'heavy', 'light']:
            if col_name in df.columns:
                expressions.append(
                    pl.col(col_name).fill_null('').str.split(by='|').list.eval(pl.element().str.strip_chars()).list.eval(pl.element().filter(pl.element() != ""))
                )
        
        mapping_list = df.with_columns(expressions).to_dicts()
        return {row['pdb_id']: {key: row.get(key, []) for key in ['antigen', 'heavy', 'light']} for row in mapping_list if row.get('pdb_id')}
    except Exception as e:
        warnings.warn(f"Could not load or process chain mapping file '{mapping_file_path}': {e}", UserWarning)
        return {}

def load_structure(path: Path, structure_id: str) -> BioStructure | None:
    """Loads a structure file, trying PDB and then MMCIF format."""
    try:
        return PDBParser(QUIET=True).get_structure(structure_id, str(path))
    except Exception:
        try:
            return MMCIFParser(QUIET=True).get_structure(structure_id, str(path))
        except Exception as e_cif:
            warnings.warn(f"BioPython failed to parse '{path}' with both PDB and MMCIF parsers: {e_cif}", UserWarning)
            return None

def _three_to_one(res3: Iterable[str]) -> str:
    """Converts a list of 3-letter amino acid codes to a 1-letter sequence."""
    return "".join(protein_letters_3to1.get(r.upper(), "X") for r in res3)

# ────────────────────────────────────────────────────────────────────────────────
AA_ALT = {
    "MSE": "M",  # selenomethionine
    "SEC": "U",  # selenocysteine
    "PYL": "O",  # pyrrolysine
    "HYP": "P",  # hydroxy‑proline
    # add more exotic codes if you need them
}

def _three_to_one(res3: Iterable[str]) -> str:
    """Convert 3‑letter codes to a 1‑letter sequence; unknown → X."""
    return "".join(
        protein_letters_3to1.get(r, AA_ALT.get(r, "X")) for r in (s.upper() for s in res3)
    )

# ────────────────────────────────────────────────────────────────────────────────
def write_valid_mmcif(structure, output_path, *, entry_id: str):
    """
    BioPython‑only mmCIF writer with explicit `stop_` and `#` separators
    so that picky parsers (e.g. PyMOL) accept the file.
    """
    # 1 ▶ write coordinates once
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(output_path), preserve_atom_numbering=True)

    # 2 ▶ collect polymer chains + sequences
    chains = []
    model  = next(structure.get_models())
    for ch in model:
        residues = [r for r in ch if r.id[0] == " "]
        if not residues:
            continue
        print(residues)
        try:
            seq = "".join(seq1(r.resname) for r in residues)
        except KeyError:
            seq = "?"
        chains.append((ch.id, seq))

    # 3 ▶ append metadata (each loop closed by `stop_` + `#`)
    today = date.today().isoformat()
    with open(output_path, "a") as fh:
        fh.write(f"""
#
_entry.id {entry_id}
_pdbx_database_status.deposit_date {today}
_pdbx_database_status.rel_release_date {today}
#
loop_
_entity.id
_entity.type
_entity.pdbx_description
""")

        for cid, _ in chains:
            fh.write(f"{cid} polymer 'Chain {cid}'\n")
        fh.write("stop_\n#\n")

        fh.write("""loop_
_struct_asym.id
_struct_asym.entity_id
""")
        for cid, _ in chains:
            fh.write(f"{cid} {cid}\n")
        fh.write("stop_\n#\n")

        fh.write("""loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code
""")
        for cid, seq in chains:
            fh.write(f"{cid} 'polypeptide(L)'\n;\n{seq}\n;\n")
        fh.write("stop_\n#\n")

    print(f"✅  CIF written (PyMOL‑compatible) → {output_path}")

def transform_chains_for_template(source_structure: BioStructure, chains_to_keep: list[str]) -> BioStructure:
    """Creates a new Bio.PDB.Structure object containing only specified chains."""
    new_structure = BioStructure(source_structure.id)
    new_model = BioModel(0)
    new_structure.add(new_model)
    source_model = next(source_structure.get_models())
    
    for chain_id in chains_to_keep:
        if source_model.has_id(chain_id):
            new_model.add(source_model[chain_id].copy())
    return new_structure

def create_residue_map(structure: BioStructure, chains_in_fasta: dict, chain_mapping_for_pdb: dict) -> dict:
    """Creates a mapping from (chain_id, res_id) to 1-based FASTA sequence index via alignment."""
    res_map = {}
    aligner = Align.PairwiseAligner(mode='global')
    pdb_to_fasta_key = {p_chain: key for key, p_chains in chain_mapping_for_pdb.items() for p_chain in p_chains}
    
    for chain in next(structure.get_models()):
        fasta_key = pdb_to_fasta_key.get(chain.id)
        if not fasta_key or fasta_key not in chains_in_fasta: continue

        pdb_residues = [res for res in chain if res.id[0] == ' ']
        try:
            pdb_seq = "".join(seq1(res.get_resname()) for res in pdb_residues)
        except KeyError:
            continue # Skip chains with non-standard residues that seq1 fails on
            
        if not pdb_seq: continue

        try:
            alignment = next(aligner.align(pdb_seq, chains_in_fasta[fasta_key]))
        except StopIteration:
            continue # Alignment failed

        pdb_res_idx, fasta_seq_idx = 0, 0
        aligned_pdb_seq, aligned_fasta_seq = alignment[0], alignment[1]

        for pdb_char, fasta_char in zip(aligned_pdb_seq, aligned_fasta_seq):
            if pdb_res_idx < len(pdb_residues):
                if pdb_char != '-':
                    if fasta_char != '-':
                        res_id_tuple = pdb_residues[pdb_res_idx].id
                        res_map[(chain.id, res_id_tuple[1])] = fasta_seq_idx + 1
                    pdb_res_idx += 1
            if fasta_char != '-':
                fasta_seq_idx += 1
    
    return res_map

def parse_fasta_for_chains(fasta_p: Path) -> dict:
    """Parses a FASTA file for antigen, heavy, and light chains based on header descriptions."""
    chains = {}
    for record in SeqIO.parse(fasta_p, "fasta"):
        h_lower = record.description.lower()
        if "antigen" in h_lower: chains["antigen"] = str(record.seq)
        elif "antibody h" in h_lower: chains["heavy"] = str(record.seq)
        elif "antibody l" in h_lower: chains["light"] = str(record.seq)
    return chains

def calculate_contacts(structure: BioStructure, heavy_chain_ids: list, antigen_chain_ids: list, antibody_chain_ids: list, residue_map: dict) -> dict | None:
    """Calculates contacts between antigen and antibody chains."""
    model = next(structure.get_models())
    antigen_atoms = [atom for cid in antigen_chain_ids if model.has_id(cid) for atom in model[cid].get_atoms() if atom.get_parent().id[0] == ' ']
    antibody_atoms = [atom for cid in antibody_chain_ids if model.has_id(cid) for atom in model[cid].get_atoms() if atom.get_parent().id[0] == ' ']
    if not antigen_atoms or not antibody_atoms: return None

    ns = NeighborSearch(antigen_atoms)
    contact_residues = set()
    dist_thresh = 4.5

    for ab_atom in antibody_atoms:
        for ag_atom in ns.search(ab_atom.coord, dist_thresh, 'A'):
            ag_res = ag_atom.get_parent()
            fasta_idx = residue_map.get((ag_res.get_parent().id, ag_res.id[1]))
            if fasta_idx: contact_residues.add((ag_res.get_parent().id, int(fasta_idx)))
    
    if not contact_residues: return None
    
    return { "pocket": { "binder": str(heavy_chain_ids[0]) if heavy_chain_ids else 'H', "contacts": ContactList([ContactPair(list(c)) for c in sorted(contact_residues)]), "max_distance": dist_thresh } }


# ==============================================================================
# --- Main Worker Function ---
# ==============================================================================

def process_file_for_all_scenarios(fasta_path_str: str, scenarios: list[Scenario], base_args: argparse.Namespace, chain_mapping: dict):
    """Processes a single FASTA file for all specified scenarios."""
    # Setup YAML representers for this worker process
    def represent_inline_list(dumper, data): return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
    yaml.add_representer(ContactList, represent_inline_list, Dumper=yaml.SafeDumper)
    yaml.add_representer(InlineList, represent_inline_list, Dumper=yaml.SafeDumper)
    def represent_contact_pair(dumper, data): return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True) if len(data) == 2 else dumper.represent_list(data)
    yaml.add_representer(ContactPair, represent_contact_pair, Dumper=yaml.SafeDumper)

    try:
        fasta_path = Path(fasta_path_str)
        pdb_id = fasta_path.stem.removesuffix('complex')
        print(f"➡️  Processing: {fasta_path.name}")

        chains = parse_fasta_for_chains(fasta_path)
        if not ("antigen" in chains and "heavy" in chains):
            print(f"    ⚠️  Skipping {fasta_path.name}: FASTA must contain 'antigen' and 'heavy' chains.")
            return

        mapping = chain_mapping.get(pdb_id, {})
        antigen_chains = [c for c in mapping.get('antigen', []) if is_valid_chain(c)]
        heavy_chains = [c for c in mapping.get('heavy', []) if is_valid_chain(c)]
        light_chains = [c for c in mapping.get('light', []) if is_valid_chain(c)]
        
        structure, constraint, cif_path_for_template = None, None, None
        needs_structure = any(s.use_template or s.use_pocket for s in scenarios)

        if needs_structure and base_args.pdb_dir:
            pdb_path = base_args.pdb_dir / f"{pdb_id}_filtered.pdb"
            if pdb_path.is_file():
                structure = load_structure(pdb_path, pdb_id)
        
        if structure:
            residue_map = create_residue_map(structure, chains, mapping)
            if any(s.use_template for s in scenarios) and antigen_chains and (heavy_chains or light_chains):
                cif_path_for_template = base_args.shared_cif_dir / f"{pdb_id}.cif"
                if not cif_path_for_template.exists():
                    template_struct = transform_chains_for_template(structure, antigen_chains + heavy_chains + light_chains)
                    write_valid_mmcif(template_struct, cif_path_for_template, entry_id=pdb_id)
            if any(s.use_pocket for s in scenarios) and residue_map:
                constraint = calculate_contacts(structure, heavy_chains, antigen_chains, heavy_chains + light_chains, residue_map)
        
        for scenario in scenarios:
            yaml_data = {"version": 1, "sequences": []}
            if chains.get("antigen"): yaml_data["sequences"].extend([{"protein": {"id": str(c), "sequence": chains["antigen"]}} for c in antigen_chains])
            if chains.get("heavy"): yaml_data["sequences"].extend([{"protein": {"id": str(c), "sequence": chains["heavy"]}} for c in heavy_chains])
            if chains.get("light"): yaml_data["sequences"].extend([{"protein": {"id": str(c), "sequence": chains["light"]}} for c in light_chains])

            if not yaml_data["sequences"]: continue

            if not scenario.use_msa:
                for s in yaml_data["sequences"]: s["protein"]["msa"] = "empty"
            if scenario.use_template and cif_path_for_template and cif_path_for_template.exists():
                template_chains = antigen_chains + heavy_chains + light_chains
                yaml_data["templates"] = [{"cif": str(cif_path_for_template.resolve()), "chain_id": InlineList(template_chains), "template_id": InlineList(template_chains)}]
            if scenario.use_pocket and constraint:
                yaml_data["constraints"] = [constraint]

            output_yaml_path = scenario.output_dir / f"{fasta_path.stem}.yaml"
            with open(output_yaml_path, 'w') as f:
                yaml.dump(yaml_data, f, Dumper=yaml.SafeDumper, sort_keys=False, indent=2, default_flow_style=False, width=1000)
        
    except Exception as e:
        warnings.warn(f"FATAL ERROR while processing {fasta_path_str}: {e}", UserWarning)
        import traceback
        traceback.print_exc()

# ==============================================================================
# --- Main Execution Logic ---
# ==============================================================================

def main():
    """Main function to parse arguments and orchestrate scenario generation."""
    parser = argparse.ArgumentParser(description="Generate YAML configuration files for antibody-antigen modeling.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # --- All 8 scenarios are defined here ---
    all_scenarios = {
        'no_template_no_msa': Scenario('no_template_no_msa', use_msa=False, use_template=False, use_pocket=False),
        'no_template_with_msa': Scenario('no_template_with_msa', use_msa=True, use_template=False, use_pocket=False),
        'no_template_no_msa_with_pocket': Scenario('no_template_no_msa_with_pocket', use_msa=False, use_template=False, use_pocket=True),
        'no_template_with_msa_with_pocket': Scenario('no_template_with_msa_with_pocket', use_msa=True, use_template=False, use_pocket=True),
        'template_no_msa_no_pocket': Scenario('template_no_msa_no_pocket', use_msa=False, use_template=True, use_pocket=False),
        'template_no_msa_with_pocket': Scenario('template_no_msa_with_pocket', use_msa=False, use_template=True, use_pocket=True),
        'template_with_msa_no_pocket': Scenario('template_with_msa_no_pocket', use_msa=True, use_template=True, use_pocket=False),
        'template_with_msa_with_pocket': Scenario('template_with_msa_with_pocket', use_msa=True, use_template=True, use_pocket=True),
    }

    parser.add_argument("--fasta_list", required=True, help="Text file listing paths to FASTA files.")
    parser.add_argument("--output_dir", type=Path, default="./benchmark_runs", help="Base directory to save all generated scenario files.")
    parser.add_argument("--pdb_dir", type=Path, help="Directory with filtered PDB files to use as templates.")
    parser.add_argument("--chain_mapping", type=Path, help="Path to chain mapping CSV file.")
    parser.add_argument("--scenarios", nargs='+', choices=all_scenarios.keys(), default=list(all_scenarios.keys()), help="Which modeling scenarios to generate.")
    parser.add_argument("--workers", type=int, default=os.cpu_count(), help="Number of worker processes to use.")
    args = parser.parse_args()

    # --- Setup output directories ---
    selected_scenarios = [all_scenarios[name] for name in args.scenarios]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.shared_cif_dir = args.output_dir / "shared_cifs"
    args.shared_cif_dir.mkdir(exist_ok=True)

    final_scenarios_to_run = []
    for scenario in selected_scenarios:
        if (scenario.use_template or scenario.use_pocket) and not args.pdb_dir:
            print(f"⏭️  Skipping scenario '{scenario.name}' because --pdb_dir is not provided.")
            continue
        scenario.output_dir = args.output_dir / "yamls" / scenario.name
        scenario.output_dir.mkdir(parents=True, exist_ok=True)
        final_scenarios_to_run.append(scenario)

    if not final_scenarios_to_run:
        print("No scenarios to run. Exiting."); return

    # --- Load input data ---
    with open(args.fasta_list) as f:
        fasta_paths = [line.strip() for line in f if line.strip()]
    
    chain_mapping = load_chain_mapping(args.chain_mapping) if args.chain_mapping and Path(args.chain_mapping).is_file() else {}
    print(f"Loaded {len(fasta_paths)} FASTA files and {len(chain_mapping)} chain mappings.")
    print(f"Generating YAMLs for {len(final_scenarios_to_run)} scenarios...")
    
    # --- Create a partial function with fixed arguments for mapping ---
    process_func = partial(process_file_for_all_scenarios, scenarios=final_scenarios_to_run, base_args=args, chain_mapping=chain_mapping)

    # --- Run processing ---
    # Multiprocessing is disabled per user request for debugging.
    # To re-enable, change the condition `if False:` to `if args.workers > 1 and len(fasta_paths) > 1:`
    if False:
        print(f"\n[INFO] Multiprocessing ENABLED. Processing files with {args.workers} workers.")
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            # Using list() to ensure all tasks complete before the program exits
            list(executor.map(process_func, fasta_paths))
    else:
        from tqdm import tqdm
        print("\n[INFO] Multiprocessing is DISABLED. Processing files sequentially.")
        for fasta_path in tqdm(fasta_paths, desc="Processing files"):
            process_func(fasta_path)

    print("\n✅  Scenario generation complete.")

if __name__ == "__main__":
    # Set start method to 'spawn' for compatibility on Linux/macOS
    try:
        if os.name != 'nt': multiprocessing.set_start_method('spawn', force=True)
    except (RuntimeError, ValueError): pass
    
    # Allow Polars to be used in forked processes
    os.environ["POLARS_ALLOW_FORKING"] = "1"
    
    main()