using QMPS

# Creating Cluster-state QMPS circuit.
nBitT1 = 4
vBit1 = 1
rBit1 = 1
mpsCS = MPSC("CS", nBitT1, vBit1, rBit1)

# Creating Differentiable QMPS circuit. 
nBitT2 = 4
vBit2 = 2
rBit2 = 1
depth2 = 2
mpsDC = MPSC(("DC", depth2), nBitT2, vBit2, rBit2)

# Calculating quantum differentiations.
# using Pkg; Pkg.add("Yao") # Un-comment this line if you haven't installed package Yao.
using Yao
nBitT0 = 6
depth0 = 4
# Differentiable circuit
c0 = DCbuilder(nBitT0, depth0).body |> markDiff
# QMPS circuit
c1 = deepcopy(mpsDC.circuit)
# MPS extended circuit
c2 = deepcopy(mpsDC.cExtend)
# Witness operators
op0 = put(nBitT0, 1=>X)
op1 = put(nqubits(c1), 1=>X)
op2 = put(nqubits(c2), 1=>X)
# Differentiable gates
dGates0 = collect_blocks(QMPS.QDiff, c0)
dGates1 = collect_blocks(QMPS.QDiff, c1)
dGates2 = collect_blocks(QMPS.QDiff, c2)
# Default registers
reg0 = zero_state(nqubits(c0))
reg1 = zero_state(nqubits(c1), nbatch=50000)
reg2 = zero_state(nqubits(c2)) 
# Quantum differentiations
grads = getQdiff!.( ()->(reg0 |> c0), dGates0, Ref(op0) )
grads = getQdiff!.( ()->(reg1 |> c1), dGates1, Ref(op1) )
grads = getQdiff!.( ()->(reg2 |> c2), dGates2, Ref(op2) )