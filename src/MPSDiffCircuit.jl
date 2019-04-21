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
        if (nBitA - vBit -rBit) % rBit != 0 
            println("Error: nBlock is not integer!")
            return 0
        end
        nBlock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        #println("M1\n")
        MPS = MPSbuilder(nBitA, nBit, nBlock, MPSblock, vBit, rBit)
        #println("M2\n")
        circuit = MPS[1]
        #println("M3\n")
        cExtend = MPS[2]
        diffs = collect(circuit, AbstractDiff)
        new(circuit, circuit.blocks, cExtend, cExtend.blocks, diffs, nBit, nBlock)
    end

end

struct DCbuilder
    block      # The Block for onr single depth.
    Cblocks    # The whole Blocks for n depths.
    body       # The differentiable circuit for n depths. 
    head       # The "head" of the circuit(for direct input qubits).
    tail       # The "tail" of the circuit(adding 2 layers of rotaional gates).
    function DCbuilder(nBit::Int64, ndepth::Int64)
        body = random_diff_circuit(nBit, ndepth, pair_ring(nBit))
        cTemp = random_diff_circuit(nBit, 2, pair_ring(nBit))
        block = chain(nBit, cTemp.blocks[3:4])
        Cblocks = chain(nBit, repeat([block],ndepth) )
        head = chain(nBit, cTemp.blocks[1:2])
        tail = cTemp.blocks[end]
        new(block, Cblocks, body, head, tail)
    end
end

# Missing cExtend!!
function MPSbuilder(nBitA::Int64, nBit::Int64, nBlock::Int64, 
                   blockT::Tuple{String, Int64}, vBit::Int64, rBit::Int64) 
    if blockT[1] == "DC" && rBit == 1 # differentiable circuit.
        depth = blockT[2]
        # println("s1\n")
        swapV1 = chain(nBit, [put(nBit, (i,i+1)=>SWAP ) for i=1:   nBit-1])
        swap1V = chain(nBit, [put(nBit, (i+1,i)=>SWAP ) for i=nBit-1:-1:1])
        # println("s2\n")
        cBlockHead = chain(nBit, DCbuilder(nBit, depth).block, swapV1) |> autodiff(:QC)
        dispatch!(cBlockHead, :random)
        # println("s3\n")
        cBlock =  chain(nBit, vcat([swap1V], cBlockHead.blocks))
        # println("s4\n")
        cBlocks = []
        cEBlocks = []        
        for i=1:nBlock
            push!(cBlocks, cBlock)
            # println("s4_1\n")
            # dispatch!(cBlocks[i], :random)
            # println("s4_2\n")
            push!(cEBlocks, put(nBitA, Tuple((nBitA-i-vBit+1):(nBitA-i+1),)=>cBlockHead))
            # println("s4_3\n")
            # dispatch!(cEBlocks[i], parameters(cBlocks[i]))
        end
        cBlocks[1] = cBlockHead
        circuit = chain(nBit, cBlocks)
        cExtend = chain(nBitA, cEBlocks)
    end

    if blockT[1] == "DC" && rBit > 1 # differentiable circuit.
        depth = blockT[2]
        # println("s1\n")
        swapV1 = chain(nBit, [put(nBit, (inBit,inBit+1)=>SWAP ) for irBit=rBit:-1:1 for inBit=irBit  :(vBit+irBit-1)])
        swap1V = chain(nBit, [put(nBit, (inBit+1,inBit)=>SWAP ) for irBit=1:   rBit for inBit=(vBit+irBit-1):-1:irBit])
        # println("s2\n")
        cBlockHead = chain(nBit, DCbuilder(nBit, depth).block, swapV1) |> autodiff(:QC)
        dispatch!(cBlockHead, :random)
        # println("s3\n")
        cBlock =  chain(nBit, vcat([swap1V], cBlockHead.blocks))
        # println("s4\n")
        cBlocks = []
        cEBlocks = []        
        for i=1:nBlock
            push!(cBlocks, cBlock)
            # println("s4_1\n")
            push!(cEBlocks, put(nBitA, Tuple((nBitA-i*rBit-vBit+1):(nBitA-(i-1)*rBit),)=>cBlockHead))
            # println("s4_2\n")
        end
        cBlocks[1] = cBlockHead
        circuit = chain(nBit, cBlocks)
        cExtend = chain(nBitA, cEBlocks)
    end

    (circuit, cExtend)
end

function MPSbuilder(nBitA::Int64, nBit::Int64, nBlock::Int64, 
                    blockT::String, vBit::Int64, rBit::Int64)

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
        circuit = chain(nBit,cBlocks)  
        cExtend = chain(nBitA, vcat([repeat(nBitA,H,Tuple(nBitA:-1:1,))], [control(nBitA, i-1, i=>Z) for i=nBitA:-1:2])) 
    end
    (circuit, cExtend)
end

end
