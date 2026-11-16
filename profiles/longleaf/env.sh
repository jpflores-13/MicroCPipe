## Site setup for UNC Longleaf, sourced by MicroCPipe.sh and unlock.sh.
## Needs python3 >= 3.11 and conda >= 24.7.1 (Snakemake 8.27): install your own Miniforge
## once (~0.5 GB) and set MICROC_CONDA_ROOT if it is not at $HOME/miniforge3.
source "${MICROC_CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate base
