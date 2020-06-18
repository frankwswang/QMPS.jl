# MPSCircuit

[![Build Status](https://travis-ci.com/frankwswang/MPSCircuit.jl.svg?branch=master)](https://travis-ci.com/frankwswang/MPSCircuit.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/frankwswang/MPSCircuit.jl?svg=true)](https://ci.appveyor.com/project/frankwswang/MPSCircuit-jl)
[![Coverage](https://codecov.io/gh/frankwswang/MPSCircuit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/frankwswang/MPSCircuit.jl)

A quick realization of MPS(Matrix Product State)-qubit circuit. 

## Supported types of MPS
- Cluster state
- Differentiable circuit constructed state

## Fields of the struct MPSC
Fields | Meanings
------------ | -------------
__circuit__|MPS circuit.
__mpsBlocks__|Array of all the MPS blocks in MPS circuit.
__cExtend__|The MPS circuit extended back to where it doesn't reuse any qubit.
__cEBlocks__|Array of all the MPS blocks in the Extended circuit.
__dGates__|Differentiable gates of the MPS circuit if applicable.
__nBit__|Number of lines(bits) of the MPS circuit. 
__nBlock__|Number of blocks in the MPS ciruict.

## Setup Guide
### Julia Environment
* [__Julia 1.3-1.4__](https://julialang.org)

### Installation
Please type `]` in Julia REPL to enter [`Pkg` mode](https://julialang.github.io/Pkg.jl/v1.0/index.html), then type:
```
pkg> add https://github.com/frankwswang/MPSCircuit.jl.git
``` 
__ATTENTION:__ This packge is dependent on package [__Yao__](https://github.com/QuantumBFS/Yao.jl) and currently compatiple version is __Yao 0.4.1__. For the future development, you need to check its compatibility if you want to use it with a higher version of __Yao__. 

## Reference
* Mitarai, K., Negoro, M., Kitagawa, M., & Fujii, K. (2018). Quantum circuit learning. Physical Review A, 98(3), 032309. ([PDF](https://arxiv.org/pdf/1803.00745.pdf))

* Liu, J. G., Zhang, Y. H., Wan, Y., & Wang, L. (2019). Variational Quantum Eigensolver with Fewer Qubits. arXiv preprint, [arXiv:1902.02663](https://arxiv.org/abs/1902.02663). ([PDF](https://arxiv.org/pdf/1902.02663.pdf))

## License
MPSCircuit.jl is released under Apache License 2.0.