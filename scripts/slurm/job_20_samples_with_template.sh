#!/bin/bash

#SBATCH -J boltz_fv_vhhs
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

#SBATCH --account SORMANNI-SL2-GPU

# Set up all cache directories in fast local storage
USER_NAME="cp864"

SCRATCH_DIR="/spinning1/scratch/${USER_NAME}/"
BOLTZ_DATA_CACHE="/spinning1/${USER_NAME}/bin/boltz_cache" # Persistent cache for Boltz's own data (CCD, etc.)
PERSISTENT_TORCH_CACHE="/spinning1/${USER_NAME}/bin/torch_cache" # Persistent cache for compiled kernel:

BOLTZ_CACHE="${SCRATCH_DIR}/boltz_cache"
BOLTZ1_CACHE="${SCRATCH_DIR}/boltz1_cache"

TORCH_CACHE_DIR="${SCRATCH_DIR}/torch_cache"
INDUCTOR_CACHE="${SCRATCH_DIR}/torch_cache/inductor"
TRITON_CACHE="${SCRATCH_DIR}/torch_cache/triton"
TRITON_CACHE_DIR="${SCRATCH_DIR}/torch_cache/triton"
CUEQ_CACHE="${SCRATCH_DIR}/torch_cache/cueq"
TORCH_COMPILE_CACHE_DIR="${SCRATCH_DIR}/torch_cache/compile/"

source /spinning1/scratch/cp864/miniforge3/etc/profile.d/conda.sh
conda activate incito-pipeline

#### FOR EXAMPLE
time boltz predict /home/${USER_NAME}/repos/incito-pipeline/incito_pipeline/datasets/out/sample_run_20/sample_yaml \
--out_dir /home/${USER_NAME}/repos/incito-pipeline/incito_pipeline/datasets/out/sample_run_20/boltz_output \
--cache "$BOLTZ_CACHE" \
 --use_msa_server --override \
 --diffusion_samples 100 \
 --max_parallel_samples 20 --preprocessing-threads 4 \
| while IFS= read -r line; do
printf "%(%Y-%m-%d %H:%M:%S)T %s\n" -1 "$line"
done | tee -a vanilla_msa_${MODEL:-boltz2}_speed.log.txt

