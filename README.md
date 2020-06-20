# QMPS

[![Build Status](https://travis-ci.com/frankwswang/QMPS.jl.svg?branch=master)](https://travis-ci.com/frankwswang/QMPS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/frankwswang/QMPS.jl?svg=true)](https://ci.appveyor.com/project/frankwswang/QMPS-jl)
[![Coverage](https://codecov.io/gh/frankwswang/QMPS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/frankwswang/QMPS.jl)

A quick realization of Quantum circuit (QMPS) of Matrix Product States compatible with a quantum simulation framework called [__Yao__](https://github.com/QuantumBFS/Yao.jl).

## Supported types of MPS

- Cluster state
- Differentiable circuit constructed state (Support quantum differentiations.)

## Main functions

### `MPSC`

Generate the structure of elements related to a QMPS circuit.

### `MPSpar`

Construct parameters that MPSC needs.

### `MPSbuilder`

Function for creating different types of MPS circuits.

### `DCbuilder`

Generate the structure of elements may needed for a Quantum differentiable circuit.

### `MPSDCpar`

Get the circuit parameters of a differentiable QMPS circuit (QMPS-DC) or of a QMPS-DC extended circuit.

### `markDiff`

Return the differentiable gate(s) `QDiff{GT, N}` from a block or a block tree such as `ChainBlock`.

### `getQdiff`

[Quantum Operator differentiation.](#jump)

### `getNdiff`

Numerical Operator differentiation.

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

## Diffrentiable QMPS circuit
### Using `expect` from 
```


```

## Setup Guide
### Julia Environment
* [__Julia 1.3-1.4__](https://julialang.org)

### Installation
Please type `]` in Julia REPL to enter [`Pkg` mode](https://julialang.github.io/Pkg.jl/v1.0/index.html), then type:
```
pkg> add https://github.com/frankwswang/MPSCircuit.jl
``` 
__ATTENTION:__ This packge is dependent on package [__Yao__](https://github.com/QuantumBFS/Yao.jl) and currently compatiple version is __Yao 0.4.1__. For the future development, you need to check its compatibility if you want to use it with a higher version of __Yao__. 

## Reference
<span id="jump">
</span>
* Mitarai, K., Negoro, M., Kitagawa, M., & Fujii, K. (2018). Quantum circuit learning. Physical Review A, 98(3), 032309. ([PDF](https://arxiv.org/pdf/1803.00745.pdf))


* Liu, J. G., Zhang, Y. H., Wan, Y., & Wang, L. (2019). Variational Quantum Eigensolver with Fewer Qubits. arXiv preprint, [arXiv:1902.02663](https://arxiv.org/abs/1902.02663). ([PDF](https://arxiv.org/pdf/1902.02663.pdf))

## License
__QMPS__ is released under Apache License 2.0.