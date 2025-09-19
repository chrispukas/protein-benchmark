#!/usr/bin/env python3
import os

from datetime import datetime

import incito_pipeline.util.parse.fasta
import rcsbapi.sequence as sequence
import incito_pipeline.util.data_utility as data_utility

from importlib.resources import files
from incito_pipeline.util import data_utility

DATASET_PATH = files("incito_pipeline.datasets")
CUTOFF_DATE = datetime(2023, 1, 6)

file_names_VHHs = data_utility.get_file_names(
    os.path.join(DATASET_PATH, "AF3_independent_test_VHHs_GroundTruth")
)
file_names_Fv = data_utility.get_file_names(
    os.path.join(DATASET_PATH, "AF3_independent_test_Fv_GroundTruth")
)

merged_paths = data_utility.merge_lists(file_names_VHHs, file_names_Fv)

"""
data_utility.build_fastas_from_pdb(os.path.join(DATASET_PATH, "AF3_independent_test_VHHs_GroundTruth"),
             os.path.join(DATASET_PATH, "out/FASTAS_VHHs"), single_output=False, cutoff_date=CUTOFF_DATE)
data_utility.build_fastas_from_pdb(os.path.join(DATASET_PATH, "AF3_independent_test_Fv_GroundTruth"),
             os.path.join(DATASET_PATH, "out/FASTAS_Fv"), single_output=False, cutoff_date=CUTOFF_DATE)
"""

#data_utility.build_fastas_from_pdb(os.path.join(DATASET_PATH, "AF3_independent_test_VHHs_GroundTruth"),
#             os.path.join(DATASET_PATH, "out/FASTAS_VHHs"), single_output=True, cutoff_date=CUTOFF_DATE)
#data_utility.build_fastas_from_pdb(os.path.join(DATASET_PATH, "AF3_independent_test_Fv_GroundTruth"),
#             os.path.join(DATASET_PATH, "out/FASTAS_Fv"), single_output=True, cutoff_date=CUTOFF_DATE)

incito_pipeline.util.parse.fasta.build_fastas_from_pdb(os.path.join(DATASET_PATH, "AF3_independent_test_VHHs_GroundTruth"),
            os.path.join(DATASET_PATH, "out/run_full/vhhs/fastas_vhhs"), 
            single_output=True, cutoff_date=CUTOFF_DATE)

