import os
import json

"""
Utility functions for data handling in the Incito pipeline.
"""

def stratize_data(dataset, acceptable_cutoff=0.23, medium_cutoff=0.49, high_cutoff=0.8):
    low_cutoff_text = f"DockQ<{acceptable_cutoff}"
    acceptable_cutoff_text = f"{acceptable_cutoff}<=DockQ<={medium_cutoff}"
    medium_cutoff_text = f"{medium_cutoff}<DockQ<{high_cutoff}"
    high_cutoff_text = f"DockQ>={high_cutoff}"

    dict_data = {
        low_cutoff_text: [],
        acceptable_cutoff_text: [],
        medium_cutoff_text: [],
        high_cutoff_text: [],
    }
            
    for i, data in enumerate(dataset):
        if float(data) >= float(high_cutoff):
            dict_data[high_cutoff_text].append(data)
        elif float(data) >= float(medium_cutoff):
            dict_data[medium_cutoff_text].append(data)
        elif float(data) >= float(acceptable_cutoff):
            dict_data[acceptable_cutoff_text].append(data)
        else:
            dict_data[low_cutoff_text].append(data)

    return dict_data

def stratize_top_n_dockq_scores(dataset, n=[1, 5, 10, 20, 50, 100]):
    sorted_scores = sorted(dataset, reverse=True)

    dict_data = {}
    for top_n in n:
        dict_data[f"top_{top_n}"] = sorted_scores[:top_n]

    return dict_data



def get_seperation_between_datasets():
    pass

def pair_up_seperations(sep_1, sep_2):
    pass


"""Utility functions to handle PDB to PDB conversion and file management."""

def get_pdb_id(file_path):
    return os.path.basename(file_path).split('_')[0]
def get_full_pdb_id(file_path):
    return os.path.basename(file_path).split('.')[0]

"""
General utility functions
"""

def merge_lists(list_1, list_2):
    return list(set(list_1) | set(list_2))

def get_file_names(dir=None, filter=None):

    if dir is None or not dir.strip():
        print("Error: No directory path provided.")
        return []
    
    if not os.path.exists(dir):
        print(f"Error: File path does not exist: {dir}")
        return []
    
    # Filter path by files only
    dirs = os.listdir(dir)
    if filter:
        dirs = [d for d in dirs if filter in d]
    return dirs
def get_dirs(dir=None, filter=None):
    if dir is None:
        raise ValueError("No directory path provided.")
    if not os.path.exists(dir):
        raise ValueError(f"File path does not exist: {dir}")
    
    # Filter path by directories only
    dirs = os.listdir(dir)
    dirs = [d for d in dirs if os.path.isdir(os.path.join(dir, d))]

    # If a filter is provided, apply it to the directory names
    if filter:
        dirs = [d for d in dirs if filter in d]

    return dirs

# Pairs up model and native PDB files from two directories.
def pairs_pdbs(model_dir, native_dir):
    is_model_dir = os.path.isdir(model_dir)
    is_native_dir = os.path.isdir(native_dir)

    # If both paths are not directories, return a dictionary with model_dir as key and native_dir as value
    if not is_model_dir and not is_native_dir:
        return {model_dir: native_dir}

    if not is_model_dir and is_native_dir:
        files_native = set(get_file_names(native_dir))
        return {
            os.path.join(native_dir, f): model_dir
            for f in set(get_file_names(native_dir))
        }
    elif is_model_dir and not is_native_dir:
        files_model = set(get_file_names(model_dir))
        return {
            os.path.join(model_dir, f): native_dir
            for f in set(get_file_names(model_dir))
        }

    files_model = set(get_file_names(model_dir))
    files_native = set(get_file_names(native_dir))
    print(files_model)
    print(files_native)

    common = files_model & files_native  # Set intersection is O(min(m,n))

    print(common)

    return {
        os.path.join(model_dir, f): os.path.join(native_dir, f)
        for f in common  # No need to sort unless output order matters
    }

"""
Utility functions for JSON reading and writing.
"""

def read_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data