#!/bin/bash

set -e

env_name="incito_pipeline"

eval "$(conda shell.bash hook)"

# Install mamba if not present
if ! command -v mamba &>/dev/null; then
    conda install mamba -n base -c conda-forge
fi

conda env remove -n $env_name -y || true
mamba env create -n $env_name --file environment.yml
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install https://github.com/openmm/pdbfixer/archive/refs/tags/1.9.tar.gz

# Activate env
conda activate $env_name

# Set LD_LIBRARY_PATH activation/deactivation scripts
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH" > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "unset LD_LIBRARY_PATH" > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

# (Optional) pip install your package itself if not editable
pip install .
