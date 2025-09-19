#!/usr/bin/env bash

ENVIRONMENT_NAME="protein_benchmark_env"


if ! command -v conda &> /dev/null; then
    echo "Conda could not be found. Please install Conda first."
    exit 1
fi


source "$(conda info --base)/etc/profile.d/conda.sh"

# Create conda environment with required dependencies
conda env list | grep -q "^${ENVIRONMENT_NAME}\s" || conda create --name "$ENVIRONMENT_NAME" python=3.12 -y
conda activate "$ENVIRONMENT_NAME" || exit 1
echo "Environment '$ENVIRONMENT_NAME' is now active."

if [ -f "environment.yml" ]; then
    echo "Updating environment from environment.yml..."
    conda env update --file environment.yml --prune -y
fi


# Install local package
echo "Installing local package in editable mode..."
pip install -e . || { echo "Failed to install local package"; exit 1; }