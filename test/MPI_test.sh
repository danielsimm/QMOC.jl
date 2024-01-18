julia +1.9.3 --project -e 'using Pkg; pkg"instantiate"'
julia +1.9.3 --project -e 'using Pkg; pkg"precompile"'
mpiexec -n 8 julia +1.9.3 --project test/MPI_test.jl