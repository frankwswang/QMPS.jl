"""
A quick realization of MPS-qubit circuit. 

For compatibility please also install in Julia Package mode:
pkg> add Yao
pkg> add https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git 

Reference: 
    Jin-Guo Liu, Yihong Zhang, Yuan Wan and Lei Wang,
    Variational Quantum Eigensolver with Fewer Qubits(arXiv:1902.02663)
"""

module MPSCircuit
export setMPSpar
export DCbuilder
export MPSC
export MPSbuilder

using Yao, Yao.Blocks
using QuAlgorithmZoo

# Structure of basic parameters. 
struct setMPSpar
    nBitA  # Number of lines(bits) in the Extended circuit. 
    vBit   # Number of virtual bits in the MPS circuit.
    rBit   # Number of reusable bits in the MPS circuit.
    nBit   # Number of lines(bits) in the MPS circuit.
    nBlock # Number of  MPS blocks in the MPS circuit. 

    function setMPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
        if (nBitA - vBit -rBit) % rBit != 0 
            println("Error: nBlock is not integer!")
            return 0
        end
        nBlock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        new(nBitA, vBit, rBit, nBit, nBlock)
    end
end

# Structure of elements may needed for a Differentiable circuit. 
struct DCbuilder
    block   # The block for 1 depth.
    Cblocks # The whole Blocks for n depths.
    body    # The Differentiable circuit for n depths. 
    head    # The "head" of the Differentiable circuit(if directly input bits).
    tail    # The "tail" of the Differentiable circuit(2 layers of rotaional gates).

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

# Structure of elements MPSCircuit.jl provides.
struct MPSC
    circuit  # MPS circuit.
    cBlocks  # Array of all the MPS blocks in the MPS circuit.
    cExtend  # The MPS circuit extended back to where it doesn't reuse any qubit.
    cEBlocks # Array of all the MPS blocks in the Extended circuit.
    dGates   # Differentiable gates of the MPS circuit if applicable.
    nBit     # Number of lines(bits) of the MPS circuit. 
    nBlock   # Number of blocks in the MPS ciruict.

    function MPSC(MPSblock, nBitA::Int64, vBit::Int64, rBit::Int64=1)
        par2nd = setMPSpar(nBitA, vBit, rBit)
        nBlock = par2nd.nBlock
        nBit = par2nd.nBit
        #println("M1\n")
        MPS = MPSbuilder(nBitA, vBit, rBit, MPSblock)
        #println("M2\n")
        circuit = MPS[1]
        #println("M3\n")
        cExtend = MPS[2]
        dGates = collect(circuit, AbstractDiff)
        new(circuit, circuit.blocks, cExtend, cExtend.blocks, dGates, nBit, nBlock)
    end

end

# Methods for creating different types of MPS circuits.
function MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::Tuple{String, Int64}) 
    par2nd = setMPSpar(nBitA, vBit, rBit)
    nBlock = par2nd.nBlock
    nBit = par2nd.nBit

    if blockT[1] == "DC" # Differentiable circuit.
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

function MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::String)
    par2nd = setMPSpar(nBitA, vBit, rBit)
    nBlock = par2nd.nBlock
    nBit = par2nd.nBit

    if blockT == "CS" # MPS blocks for 1D cluster state.
        if vBit !=1 || rBit !=1
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
