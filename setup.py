#!/usr/bin/env python

from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(
        name="incito_pipeline",
        version="0.0.1",
        packages=find_packages(),
        package_data={"incito_pipeline": ["datasets/AF3_independent_test_Fv_GroundTruth/*.pdb",
                      "datasets/AF3_independent_test_VHHs_GroundTruth/*.pdb",
                      "datasets/pdb_fasta"]},
        description="incito_pipeline: A pipeline for benchmarking and analyzing protein docking models.",
    )
