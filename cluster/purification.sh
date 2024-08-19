#!/bin/bash -l
#SBATCH --exclusive
#SBATCH --nodes=30
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -p smp
#SBATCH -J "QMOC_YKFullPurificationDetail24"
#SBATCH --time=20:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=simm@thp.uni-koeln.de
#SBATCH -A ag-trebst
#SBATCH --output=YKFullPurificationDetail24.out
#SBATCH --error=YKFullPurificationDetail24.err

cd QMOC.jl/
julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
srun -n 30 julia -t192 --project cluster/YKFullPurification.jl