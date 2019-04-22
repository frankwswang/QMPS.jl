# MPSCircuit.jl
A quick realization of MPS(Matrix Product State)-qubit circuit. 

# Supported types of MPS
- Cluster state
- Differentiable circuit constructed state

# Elements of the struct MPSC
Elements | Meanings
------------ | -------------
__circuit__|MPS differentiable circuit.
__cBlocks__|Array of all the MPS blocks in MPS circuit.
__cExtend__|MPS circuit extended back to where it doesn't reuse any qubit.
__cEBlocks__|Array of all the MPS blocks in extended circuit.
__diffs__|Differentials of MPS circuit.
__nBit__|Number of lines(bits) of MPS circuit. 
__nBlock__|Number of blocks in MPS ciruict.

# Setup Guide
## Julia Environment
* [__Julia 1.1__](https://julialang.org)

## Installation
Please type `]` in Julia REPL to enter `Pkg` mode, then type:
```
pkg> add https://github.com/frankwswang/MPSDiffCircuit.jl.git
``` 

## Dependent Packges
Please use the same method as above to install the packages below:

[__Yao__ (v0.3.2)](https://github.com/QuantumBFS/Yao.jl)
```
pkg> add Yao
``` 

[__QuAlgorithmZoo__](https://github.com/QuantumBFS/QuAlgorithmZoo.jl)
```
pkg> add https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git
``` 


# Reference
* Variational Quantum Eigensolver with Fewer Qubits ([pdf](https://arxiv.org/pdf/1902.02663.pdf)), [arXiv:1902.02663](https://arxiv.org/abs/1902.02663), Jin-Guo Liu, Yihong Zhang, Yuan Wan and Lei Wang

# License
MPSDiffCircuit.jl is released under Apache License 2.0.