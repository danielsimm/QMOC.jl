#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=26
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -p smp
#SBATCH -J "QMOC_FullKitaev36"
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de
#SBATCH -A ag-trebst
#SBATCH --output=FullKitaev36.out
#SBATCH --error=FullKitaev36.err

cd QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 26 julia -t192 --project cluster/Full_Kitaev36.jl