# Quantum Measurement-Only Circuits (QMOC.jl)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://danielsimm.github.io/QMOC.jl/dev/)
[![Build Status](https://github.com/danielsimm/LatticeCircuits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/danielsimm/LatticeCircuits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/danielsimm/QMOC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/danielsimm/QMOC.jl)

Pure Julia package implementing large scale measurement-only quantum circuit simulations on different geometries. [https://github.com/QuantumSavory/QuantumClifford.jl](QuantumClifford.jl) used to implement Clifford algebra. Simulation parallelism implemented with MPI, currently without fallback - make sure you have a working MPI executable linked. 

Mostly undocumented and early stage of development, expect breaking changes and added features.
