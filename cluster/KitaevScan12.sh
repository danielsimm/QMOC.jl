#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=100
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH -p normal
#SBATCH -J "kitaev_span_128"
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de

module load mpi
module load lang
module load JuliaHPC
cd /scratch/hpc-prf-pm2frg/simm/QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 100 julia -t128 --project cluster/KitaevScan12.jl