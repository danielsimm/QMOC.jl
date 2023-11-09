julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
mpiexec -n 8 julia --project test/MPI_test.jl