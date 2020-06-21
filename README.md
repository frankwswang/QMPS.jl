# QMPS

[![Build Status](https://travis-ci.com/frankwswang/QMPS.jl.svg?branch=master)](https://travis-ci.com/frankwswang/QMPS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/frankwswang/QMPS.jl?svg=true)](https://ci.appveyor.com/project/frankwswang/QMPS-jl)
[![Coverage](https://codecov.io/gh/frankwswang/QMPS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/frankwswang/QMPS.jl)

A quick realization of qubit-efficient quantum circuit architecture of Matrix Product States (QMPS). This package is an extension package for a quantum simulation framework called [__Yao__](https://github.com/QuantumBFS/Yao.jl).

## Supported types of MPS

- Cluster state
- Differentiable circuit constructed state (Support quantum differentiations)

## Main functions

### For constructing QMPS circuits

*`MPSC`: Generate the structure of elements related to a QMPS circuit.

*`MPSpar`: Construct parameters that MPSC needs.

*`MPSbuilder`: Function for creating different types of MPS circuits.

### For differentiable quantum circuits

*`DCbuilder`: Generate the structure of elements may needed for a Quantum differentiable circuit.

*`MPSDCpar`: Get the circuit parameters of a differentiable QMPS circuit (QMPS-DC) or of a QMPS-DC extended circuit.

*`markDiff`: Return the differentiable gate(s) `QDiff{GT, N}` from a block or a block tree such as `ChainBlock`.

*`getQdiff`: [Quantum differentiation.](#jump)

*`getNdiff`: Numerical differentiation.

For more introductions and tutorials about the above functions please check the __examples__ directory in the repository as well as the function documentation using Julia's [__`Help` mode__](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode-1).

## Fields of struct `MPSC`

Fields | Meanings
------------ | -------------
__circuit__|QMPS circuit.
__mpsBlocks__|Array of all the MPS blocks in the QMPS circuit.
__cExtend__|The circuit QMPS circuit is extended back to that doesn't reuse any qubit.
__cEBlocks__|Array of all the MPS blocks in the Extended circuit.
__dGates__|Differentiable gates of the QMPS circuit if applicable.
__nBit__|Number of lines (bits) of the QMPS circuit.
__nBlock__|Number of blocks in the QMPS circuit.

## Setup Guide

### Environment

*[__Julia 1.3-1.4__](https://julialang.org)
*[__Yao 0.6__](https://github.com/QuantumBFS/Yao.jl)

### Installation

Please type `]` in Julia REPL to enter [__`Pkg` mode__](https://julialang.github.io/Pkg.jl/v1.0/index.html), then type:

```julia
pkg> add https://github.com/frankwswang/QMPS.jl
```

<span id="jump">
</span>

## Reference

*[Mitarai, K., Negoro, M., Kitagawa, M., & Fujii, K. (2018). Quantum circuit learning. Physical Review A, 98(3), 032309.](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.032309)

*[Liu, J. G., Zhang, Y. H., Wan, Y., & Wang, L. (2019). Variational quantum eigensolver with fewer qubits. Physical Review Research, 1(2), 023025.](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.023025)

## License

__QMPS__ is released under Apache License 2.0.
