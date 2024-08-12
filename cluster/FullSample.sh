#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=26
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -p smp
#SBATCH -J "QMOC_FullKekule72"
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de
#SBATCH -A ag-trebst
#SBATCH --output=FullKitaevNNN72.out
#SBATCH --error=FullKitaevNNN72.err

cd QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 26 julia -t192 --project cluster/Full_Kekule72.jl