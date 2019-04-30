# MPSCircuit.jl
A quick realization of MPS(Matrix Product State)-qubit circuit. 

# Supported types of MPS
- Cluster state
- Differentiable circuit constructed state

# Elements of the struct MPSC
Elements | Meanings
------------ | -------------
__circuit__|MPS circuit.
__cBlocks__|Array of all the MPS blocks in MPS circuit.
__cExtend__|The MPS circuit extended back to where it doesn't reuse any qubit.
__cEBlocks__|Array of all the MPS blocks in the Extended circuit.
__dGates__|Differentiable gates of the MPS circuit if applicable.
__nBit__|Number of lines(bits) of the MPS circuit. 
__nBlock__|Number of blocks in the MPS ciruict.

# Setup Guide
## Julia Environment
* [__Julia 1.1__](https://julialang.org)

## Installation
Please type `]` in Julia REPL to enter `Pkg` mode, then type:
```
pkg> add https://github.com/frankwswang/MPSCircuit.jl.git
``` 
__ATTENTION:__ This packge is dependent on Julia package [__Yao__ (v0.3.2)](https://github.com/QuantumBFS/Yao.jl) and is currently compatiple with __Yao__'s version __0.4.1__. For the future development, you may need to ckeck its compatibility if you find the __Yao__'s version installed is different from the confirmed compatible version. 

# Reference
* Variational Quantum Eigensolver with Fewer Qubits ([pdf](https://arxiv.org/pdf/1902.02663.pdf)), [arXiv:1902.02663](https://arxiv.org/abs/1902.02663), Jin-Guo Liu, Yihong Zhang, Yuan Wan and Lei Wang

# License
MPSCircuit.jl is released under Apache License 2.0.
