push!(LOAD_PATH,abspath("src"))

using MPSCircuit

# Creating Cluster-state MPS circuit.
nBitT = 4
vBit = 1
rBit = 1
MPS_CScircuit = MPSC("CS",nBitT,vBit,rBit).circuit
@show MPS_CScircuit

# Creating Differentiable MPS circuit. 
nBitT = 4
vBit = 2
rBit = 1
depth = 2
MPS_Dcircuit = MPSC(("DC",depth),nBitT,vBit,rBit).circuit
@show MPS_Dcircuit