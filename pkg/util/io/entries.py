from dataclasses import dataclass


@dataclass
class Template:
    cif_path: str = None
    pdb_path: str = None

    force: bool = None
    threshold: float = None

    chain_id: list[str] = None
    template_id: list[str] = None

@dataclass
class PocketConstraint:
    binder: str = None
    contacts: list[list[str]] = None
    max_distance: float = None
    force: bool = None

@dataclass
class Sequence:
    entity_type: str = None
    id: str = None
    sequence: str = None
    msa: str = None

    position_modification: int = None



# ========================
# --- Template Entries ---
# ========================

def create_template_entry(elem: Template) -> dict:
    entry = {}

    # Populate elements
    if elem.cif_path:
        entry['cif'] = str(elem.cif_path)
    elif elem.pdb_path:
        entry['pdb'] = str(elem.pdb_path)

    if elem.force:
        entry['force'] = elem.force
        entry['threshold'] = elem.threshold
    
    if elem.chain_id and len(elem.chain_id) > 0:
        entry['chain_id'] = [str(cid) for cid in elem.chain_id]
    if elem.template_id and len(elem.template_id) > 0:
        entry['template_id'] = [str(tid) for tid in elem.template_id]

    return entry

def get_template_entries(templates: list[Template]) -> list[dict]:
    entries = []
    for elem in templates:
        entry = create_template_entry(elem)
        if entry:
            entries.append(entry)
    return entries

# ==========================
# --- Constraint Entries ---
# ==========================

def create_constraint_entry(constraint):
    fmap = {PocketConstraint: create_pocket_entry}
    return fmap[type(constraint)](constraint)

def create_pocket_entry(elem: PocketConstraint) -> dict:
    entry = {}

    # Populate elements
    if elem.binder:
        entry['binder'] = str(elem.binder)
    if elem.contacts:
        entry['contacts'] = [list(contact) for contact in elem.contacts]
    if elem.max_distance is not None:
        entry['max_distance'] = elem.max_distance
    if elem.force is not None:
        entry['force'] = elem.force

    return {'pocket': entry}

def get_constraint_entries(constraints: list) -> list[dict]:
    entries = []
    for elem in constraints:
        entry = create_constraint_entry(elem)
        if entry:
            entries.append(entry)

    return entries

# ==========================
# --- Protein Entries ---
# ==========================

def create_sequence_entry(elem: Sequence) -> dict:
    if not elem.entity_type:
        return None

    # nested dict with sequence data
    seq_data = {}
    if elem.id:
        seq_data['id'] = str(elem.id)
    if elem.msa:
        seq_data['msa'] = str(elem.msa)
    if elem.sequence:
        seq_data['sequence'] = str(elem.sequence)
    if elem.position_modification:
        seq_data['modifications'] = [{'position': elem.position_modification}]

    # wrap under entity_type
    return {elem.entity_type: seq_data}

def get_sequence_entries(sequences: list[Sequence]) -> list[dict]:
    entries = []
    for elem in sequences:
        entry = create_sequence_entry(elem)
        if entry:
            entries.append(entry)

    return entries




