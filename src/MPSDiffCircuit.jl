module MPSDiffCircuit

"""
A quick realization of MPS-differentiable-qubit circuit. 

For compatibility please also install in Julia Package mode:
add Yao
add https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git 
"""
using Yao, Yao.Blocks
using QuAlgorithmZoo

export MPSDC

struct MPSDC
    circuit  # MPS differentiable circuit
    cBlocks  # Array of all the MPS blocks in MPS circuit.
    cExtend  # The circuit after extending the MPS circuit.
    cEBlocks # Array of all the MPS blocks in extended circuit.
    diffs    # differentials of MPS circuit.
    nBit     # Number of lines(bits) of MPS circuit. 
    
    function MPSDC(depth::Int64, nBitA::Int64, vBit::Int64, rBit::Int64=1)
        (nBitA - vBit -rBit) % rBit != 0 && println("Error: Nblock is not integer!")
        Nblock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        cBlock = random_diff_circuit(nBit, depth, pair_ring(nBit)) |> autodiff(:QC)
        cBlocks = []
        cEBlocks = []
        for i = 1:Nblock
            cBlocks = push!(cBlocks, cBlock)
            cEBlocks = push!(cEBlocks, put( nBitA, Tuple((nBitA-rBit*(i-1)):-1:(nBitA-rBit*i-vBit+1),)=>cBlock) )
        end
        circuit = chain(nBit,cBlocks)
        cExtend = chain(nBitA, cEBlocks) 
        circuit = dispatch!(circuit, :random)
        diffs = collect(circuit, AbstractDiff)
        new(circuit, cBlocks, cExtend, cEBlocks, diffs, nBit)
    end
end

end