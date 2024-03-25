julia --project -e 'using Pkg; pkg"instantiate"'
julia --project -e 'using Pkg; pkg"precompile"'
mpiexec -n 2 julia -t4 --project multithreaded_YK.jl