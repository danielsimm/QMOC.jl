#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -p smp
#SBATCH -J "QMOC_Kitaev88"
#SBATCH --time=100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de
#SBATCH -A ag-trebst
#SBATCH --output=Kitaev88.out
#SBATCH --error=Kitaev88.err

cd QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 50 julia -t192 --project cluster/Kitaev88.jl