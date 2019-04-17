"""
A quick realization of MPS-differentiable-qubit circuit. 

For compatibility please also install in Julia Package mode:
pkg> add Yao
pkg> add https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git 
"""

module MPSDiffCircuit
using Yao, Yao.Blocks
using QuAlgorithmZoo

export MPSDC
export MPSblocks

struct MPSDC
    circuit  # MPS differentiable circuit
    cBlocks  # Array of all the MPS blocks in MPS circuit.
    cExtend  # MPS circuit extended back to where it doesn't reuse any qubit.
    cEBlocks # Array of all the MPS blocks in extended circuit.
    diffs    # differentials of MPS circuit.
    nBit     # Number of lines(bits) of MPS circuit. 
    nBlock   # Number of blocks in MPS ciruict.
    
    function MPSDC(MPSblock, nBitA::Int64, vBit::Int64, rBit::Int64=1)
        (nBitA - vBit -rBit) % rBit != 0 && println("Error: nBlock is not integer!")
        nBlock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        cBlocks = MPSblocks(nBit, nBlock, MPSblock)
        cEBlocks = []
        for i = 1:nBlock
            cBlocks[i] |> autodiff(:QC)
            push!(cEBlocks, put( nBitA, Tuple((nBitA-rBit*i-vBit+1):(nBitA-rBit*(i-1)),)=>cBlocks[i]) )
        end
        circuit = chain(nBit,cBlocks)
        cExtend = chain(nBitA, cEBlocks) 
        circuit = dispatch!(circuit, :random)
        diffs = collect(circuit, AbstractDiff)
        new(circuit, cBlocks, cExtend, cEBlocks, diffs, nBit, nBlock)
    end

end

function MPSblocks(nBit::Int64, nBlock::Int64, blockT::Tuple{String, Int64})
    if blockT[1] == "DC" # differentiable circuit.
        depth = blockT[2]
        GATE = CNOT
        CXblocks =chain(nBit,
                    chain(nBit,[put(nBit, ( i,(i-1) )=>GATE) for i = nBit:-1:2]),
                    put(nBit, ( 1,nBit ) =>GATE)
        )
        cBlock = chain(nBit, random_diff_circuit(nBit, depth, pair_ring(nBit)), CXblocks)
        cBlocks = []        
        for i=1:nBlock
            push!(cBlocks, cBlock)
        end
    end
    cBlocks
end

function MPSblocks(nBit::Int64, nBlock::Int64, blockT::String)
    if blockT == "U1" # VQE MPS blocks for U(1) symmetry.
        cBlock = chain(nBit,
                       vcat([chain(nBit, repeat(nBit, Rz(0), (i, i-1)), put(nBit, (i, i-1)=>SWAP) ) for i=nBit:-1:2 ],
                            [chain(nBit, put(nBit, nBit=>Rz(0)), put(nBit, 1=>Rz(0)))],
                            [put(nBit, (nBit,1)=>SWAP)]
                        )
                      ) 
        cBlocks = []
        for i = 1:nBlock
            cBlockOdd = chain(nBit, put(nBit, nBit=>X), cBlock)
            cBlockEven = cBlock 
            if i%2 == 1
                push!(cBlocks, cBlockOdd)
            else
                push!(cBlocks, cBlockEven)
            end
        end
    end

    if blockT == "CS" # MPS blocks for 1D cluster state.
        if nBit != 2
            println("Error: The nBit=vBit+rBit of cluster state MPS blocks should be 2!\n")
            return 0
        end
        cBlocks = []
        push!(cBlocks, chain(nBit, repeat(nBit,H,(2,1)), control(nBit, 1, nBit=>Z)))
        for i =2:nBlock
            push!(cBlocks, chain(nBit,put(nBit, (2,1)=>SWAP), put(nBit, 1=>H), control(nBit, 1, nBit=>Z)))
        end  
    end
    cBlocks
end

end
