using Yao, Yao.Blocks
using QuAlgorithmZoo

struct MPSDC
    circuit
    cBlocks
    cExpend
    cEBlocks
    diffs
    nBit
    
    function MPSDC(depth::Int64, nBitA::Int64, vBit::Int64, rBit::Int64=1)
        (nBitA - vBit -rBit) % rBit != 0 && println("Error: Nblock is not integer!")
        Nblock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        cBlock = random_diff_circuit(nBit, depth, pair_ring(nBit)) 
        cBlocks = []
        cEBlocks = []
        for i = 1:Nblock
            cBlocks = push!(cBlocks, cBlock)
            cEBlocks = push!(cEBlocks, put( nBitA, Tuple((nBitA-rBit*(i-1)):-1:(nBitA-rBit*i-vBit+1),)=>cBlock) )
        end
        circuit = chain(nBit,cBlocks)
        cExpand = chain(nBitA, cEBlocks) |> autodiff(:QC)
        circuit = dispatch!(circuit, :random) |> autodiff(:QC)
        diffs = collect(circuit, AbstractDiff)
        new(circuit, cBlocks, cExpand, cEBlocks, diffs, nBit)
    end
end

mps = MPSDC(1,4,2)

 
