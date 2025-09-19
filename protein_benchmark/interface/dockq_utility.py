import os
import json 
import time

import numpy as np

from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

import protein_benchmark.data_utility as du
import protein_benchmark.cache as cache


from Bio.PDB import PDBParser, Chain, PDBIO


def extract_merge_chain_map_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    chain_map = []
    for model in structure:
        for chain in model:
            chain_map.append(chain.id)

    if 'L' in chain_map:
        chain_map = {
            'A': ['A'],
            'B': ['H', 'L']
        }
    else:
        chain_map = {
            'A': ['A'],
            'B': ['H']
        }
    
    return chain_map

def merge_chains(model, chain_map):
    chain_lengths = []

    # Flatten the chain_map to get all existing chains
    for new_chain_id, old_chain_ids in chain_map.items():
        target_chain = model[old_chain_ids[0]]
        model.detach_child(old_chain_ids[0])

        target_chain.id = new_chain_id
        model.add(target_chain)
        
        print(f"Mapping {new_chain_id} <- {old_chain_ids[0]}")

        for i, old_chain_id in enumerate(old_chain_ids[1:]):
            print(f"Mapping {new_chain_id} <- {old_chain_id}")

            if new_chain_id == old_chain_id:
                continue

            for j, res in enumerate(model[old_chain_id]):
                res.id = (old_chain_id, res.id[1] + i * 10, res.id[2])
                target_chain.add(res)
            model.detach_child(old_chain_id)   

        target_chain_length = len(target_chain)
        chain_lengths.append(target_chain_length)

    return model, chain_lengths

def get_dockq_score(model_pdb_path, native_pdb_path, chain_map=None, merge_chain_map=None):
    print(f"Getting DockQ Score for model: {os.path.basename(model_pdb_path)} and native: {os.path.basename(native_pdb_path)}")
    try:
        model_pdb = load_PDB(model_pdb_path)
        native_pdb = load_PDB(native_pdb_path)
    except Exception as e:
        print(f"Error loading PDB files {model_pdb_path}, {native_pdb_path}: {e}")
        return None, None

    if merge_chain_map is None:
        extracted_chain_map = extract_merge_chain_map_from_pdb(model_pdb_path)
        merge_chain_map = extracted_chain_map
        chain_map = {'A':'A', 'B':'B'}
        print("Extracted chain map from PDB:", merge_chain_map)

    model_pdb, model_chain_lengths = merge_chains(model_pdb, merge_chain_map)
    native_pdb, native_chain_lengths = merge_chains(native_pdb, merge_chain_map)

    print(f"Model residues: {model_chain_lengths}, Native residues: {native_chain_lengths}")

    try:
        full_result, dockq_result = run_on_all_native_interfaces(model_pdb, native_pdb, chain_map=chain_map)
    except Exception as e:
        io = PDBIO()

        io.set_structure(model_pdb)
        io.save("/home/cp864/repos/incito-pipeline/incito_pipeline/datasets/out/sample_run_20/failures/model_pdb.pdb")

        io.set_structure(native_pdb)
        io.save("/home/cp864/repos/incito-pipeline/incito_pipeline/datasets/out/sample_run_20/failures/native_pdb.pdb")

        print(f"Error running DockQ for {model_pdb_path}, {native_pdb_path}: {e}")
        return None, None

    return full_result, dockq_result

def get_dockq_score_from_pairs(pairs, chain_map=None, merge_chain_map=None):
    full_results = []
    dockq_scores = []
    names = []

    failed_count = 0

    for model_dir, native_dir in pairs.items():
        full_result, dockq_score = get_dockq_score(model_dir, 
                                                    native_dir, 
                                                    chain_map=chain_map, 
                                                    merge_chain_map=merge_chain_map)
        
        if full_result is None or dockq_score is None:
            failed_count += 1
            continue

        pdb_name = os.path.basename(model_dir)
        
        full_results.append(full_result)
        dockq_scores.append(dockq_score)
        names.append(pdb_name)

    return full_results, dockq_scores, failed_count, names

# FILE PATTERN MUST BE:
# <model_name>_<native_name>_<chain_map>.pdb
# PREDICTED_PDB
# └<model_name>_<native_name>_<chain_map>
#     └<model_name>_<native_name>_<chain_map>_<custom_file_naming>.pdb
#      <...>
#  <...> (continued directories) -> searches root, then using model_name, finds corresponding native PDB in this directory
#
# Native PDBs are the predicted PDBs, and the model PDBs are the ground truth PDBs.
#



def get_confidence_scores(prediction_source_dir, score_type='iptm'):
    confidence_files = du.get_file_names(prediction_source_dir, ".json")
    confidence_scores = {}

    for file in confidence_files:
        file_path = os.path.join(prediction_source_dir, file)
        data = du.read_json(file_path)
        score = data.get(score_type, 0.0)
        confidence_scores[file_path] = score

    return confidence_scores

def pick_highest_confidence(prediction_source_dir, score_type='iptm'):
    confidence_scores = get_confidence_scores(prediction_source_dir, score_type=score_type)

    if not confidence_scores:
        return None, None

    max_conf_file = max(confidence_scores, key=confidence_scores.get)
    max_conf = confidence_scores[max_conf_file]

    return max_conf, os.path.basename(max_conf_file)

def pick_lowest_confidence(prediction_source_dir):
    confidence_scores = get_confidence_scores(prediction_source_dir)
    
    if not confidence_scores:
        return None, None

    min_conf_file = min(confidence_scores, key=confidence_scores.get)
    min_conf = confidence_scores[min_conf_file]

    return min_conf, os.path.basename(min_conf_file)



def score_highest_confidence_by_dockq(path_to_gt, path_to_predicted_pdbs, path_to_boltz_output, 
                                      chain_map=None,
                                      cache_name=None,
                                      score_type='iptm'):
    gt_files = du.get_file_names(path_to_gt, filter='.pdb')

    res = []

    fail_count = 0

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        if c.is_pickle():
            res = c.load_pickle()

    if res:
        print(f"Loaded {len(res[0])} DockQ scores from cache.")
        return res
    
    for gt_file in gt_files:
        gt_name = gt_file.split('.')[0]
        gt_path = os.path.join(path_to_gt, gt_file)

        conf_dir_name = f"{gt_name}_combined"
        conf_path = os.path.join(path_to_boltz_output, conf_dir_name)

        if not os.path.exists(conf_path):
            print(f"Confidence path does not exist: {conf_path}")
            continue

        max_conf, max_conf_file = pick_highest_confidence(conf_path, score_type=score_type)
        max_conf_file_id = max_conf_file.split(gt_name)[-1].split('.')[0]
        
        pred_pdb_name = f"{gt_name}{max_conf_file_id}.pdb"
        pred_pdb_path = os.path.join(path_to_predicted_pdbs, gt_name, pred_pdb_name)

        pairs = du.pairs_pdbs(gt_path, pred_pdb_path)

        full_res, dockq, failed, names = get_dockq_score_from_pairs(pairs, chain_map=None)
        fail_count += failed

        if full_res == []:
            continue
        if len(full_res) == 0:
            full_res = full_res[0]
        if len(dockq) == 0:
            dockq = dockq[0]

        res.append((full_res, dockq, max_conf))

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        c.save_pickle(res, cache_name)


    print(f"Got {len(res)} DockQ scores, {fail_count} failed.")
    return res

def score_lowest_confidence_by_dockq(path_to_gt, path_to_predicted_pdbs, path_to_boltz_output, 
                                     chain_map=None,
                                     cache_name=None):
    gt_files = du.get_file_names(path_to_gt, filter='.pdb')

    full_results = []
    dockq_scores = []

    fail_count = 0

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        if c.is_pickle():
            full_results, dockq_scores = c.load_pickle()
    
    if full_results and dockq_scores:
        print(f"Loaded {len(full_results)} DockQ scores from cache.")
        return full_results, dockq_scores
    
    for gt_file in gt_files:
        gt_name = gt_file.split('.')[0]
        gt_path = os.path.join(path_to_gt, gt_file)

        conf_dir_name = f"{gt_name}_combined"
        conf_path = os.path.join(path_to_boltz_output, conf_dir_name)

        if not os.path.exists(conf_path):
            print(f"Confidence path does not exist: {conf_path}")
            continue

        min_conf, min_conf_file = pick_lowest_confidence(conf_path)
        min_conf_file_id = min_conf_file.split(gt_name)[-1].split('.')[0]
        
        pred_pdb_name = f"{gt_name}{min_conf_file_id}.pdb"
        pred_pdb_path = os.path.join(path_to_predicted_pdbs, gt_name, pred_pdb_name)

        pairs = du.pairs_pdbs(gt_path, pred_pdb_path)

        full_res, dockq, failed = get_dockq_score_from_pairs(pairs, chain_map=None)
        fail_count += failed

        if full_res == []:
            continue

        full_results.append(full_res)
        dockq_scores.append(dockq)

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        c.save_pickle((full_results, dockq_scores), cache_name)

    print(f"Got {len(full_results)} DockQ scores, {fail_count} failed.")
    return full_results, dockq_scores

def score_all_by_dockq_nested(path_to_gt, path_to_pred_pdbs, chain_map=None, cache_name=None, get_structure_names=True) -> dict:
    gt_files = du.get_file_names(path_to_gt, filter='.pdb')

    full_results = []
    dockq_scores = []
    structure_names = []

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        if c.is_pickle():
            full_results, dockq_scores, structure_names = c.load_pickle()

    if full_results and dockq_scores:
        print(f"Loaded {sum(len(sub) for sub in full_results)} DockQ scores from cache.")
        return full_results, dockq_scores, structure_names
    print(f"Scoring {len(gt_files)} PDB files by DockQ...")

        
    # Select the maximum DockQ score for each set of 10 predictions
    for gt_file in gt_files:
        gt_name = gt_file.split('.')[0]
        gt_path = os.path.join(path_to_gt, gt_file)
        
        pred_dir = os.path.join(path_to_pred_pdbs, gt_name)
        pairs = du.pairs_pdbs(gt_path, pred_dir)

        
        full_res, dockq, failed = get_dockq_score_from_pairs(pairs, chain_map=None)

        if full_res == []:
            continue

        full_results.append(full_res)
        dockq_scores.append(dockq)

        if not get_structure_names:
            continue
        print("Getting structure names for:", gt_name)

        structures = {}
        for i, (key, value) in enumerate(pairs.items()):
            structures[dockq[i]] = {key: value}
        structure_names.append(structures)

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        c.save_pickle((full_results, dockq_scores, structure_names), cache_name)

    print(f"Got {sum(len(sub) for sub in full_results)}")
    return full_results, dockq_scores, structure_names



def score_all_by_iptm_dockq_nested(path_to_gt, 
                                   path_to_pred_pdbs, 
                                   path_to_boltz_output, 
                                   chain_map=None, 
                                   cache_name=None,
                                   score_type="iptm") -> dict:
    gt_files = du.get_file_names(path_to_gt, filter='.pdb')

    res = []

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        if c.is_pickle():
            res = c.load_pickle()
            print(f"Loaded {sum(len(sub) for sub in res)} DockQ scores from cache.")
            return res

    for gt_file in gt_files:
        gt_name = gt_file.split('.')[0]
        gt_path = os.path.join(path_to_gt, gt_file)

        conf_dir_name = f"{gt_name}_combined"
        conf_path = os.path.join(path_to_boltz_output, conf_dir_name)

        conf_scores = get_confidence_scores(conf_path, score_type=score_type)
        print(f"Extracted {conf_scores} confidence scores for {gt_name}.")
        
        pred_dir = os.path.join(path_to_pred_pdbs, gt_name)
        pairs = du.pairs_pdbs(gt_path, pred_dir)

        full_res, dockq, failed, names = get_dockq_score_from_pairs(pairs, chain_map=None)

        if full_res == []:
            continue

        res.append((full_res, dockq, conf_scores))

    if cache_name is not None:
        c = cache.Cache(cache_dir=None, pickle_name=cache_name)
        c.save_pickle(res, cache_name)


    
    print(f"Got {sum(len(sub) for sub in res)} DockQ scores.")

    return res