#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -p smp
#SBATCH -J "QMOC_KitaevNNNedge_crit_40"
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de
#SBATCH -A ag-trebst

cd QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 2 julia -t192 --project cluster/KitaevNNN_edge_critical_40.jl