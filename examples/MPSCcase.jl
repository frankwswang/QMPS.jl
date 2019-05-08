using MPSCircuit

# Creating Cluster-state MPS circuit.
nBitT = 4
vBit = 1
rBit = 1
mpsCS = MPSC("CS",nBitT,vBit,rBit)
@show mpsCS

# Creating Differentiable MPS circuit. 
nBitT = 4
vBit = 2
rBit = 1
depth = 2
mpsDC = MPSC(("DC",depth),nBitT,vBit,rBit)
@show mpsDC