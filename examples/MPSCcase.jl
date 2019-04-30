using MPSCircuit

# Creating Cluster-state MPS circuit.
nBitT = 4
vBit = 1
rBit = 1
MPS_CScircuit = MPSC("CS",nBitT,vBit,rBit)
@show MPSGen

# Creating Differentiable MPS circuit. 
nBitT = 4
vBit = 2
rBit = 2
depth = 2
MPS_Dcircuit = MPSC(("DC",depth),nBitT,vBit,rBit)
@show MPSGen


