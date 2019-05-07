"""
Circuit builders for MPSCircuit.jl.
"""

export DCbuilder, MPSbuilder

# Structure of elements may needed for a Differentiable circuit. 
struct DCbuilder
    block   # The block for 1 depth.
    blocks  # The whole Blocks for n depths.
    Cblocks # The whole Blocks chained.
    body    # The Differentiable circuit for n depths. 
    head    # The "head" of the Differentiable circuit(if directly input bits).
    tail    # The "tail" of the Differentiable circuit(2 layers of rotaional gates).

    function DCbuilder(nBit::Int64, ndepth::Int64) #DON'T USE repeat to avoid POINTER side-effect!
        block = chain(nBit, vcat([put(nBit, i=>Rz(0)) for i=nBit:-1:1], 
                                 [put(nBit, i=>Rx(0)) for i=nBit:-1:1], 
                                 [put(nBit, i=>Rz(0)) for i=nBit:-1:1], 
                                 [control(nBit, i, (i-1)=>X) for i=nBit:-1:2], 
                                 control(nBit, 1, nBit=>X)))
        tail = chain(nBit, vcat([put(nBit, i=>Rz(0)) for i=nBit:-1:1], 
                                [put(nBit, i=>Rx(0)) for i=nBit:-1:1]))
        head = deepcopy(block[2:end])
        blocks = []
        for i=1:ndepth push!(blocks, block) end
        Cblocks = chain(nBit, blocks)
        body = chain(nBit, vcat([head], Cblocks, [tail])) 
        new(block, Cblocks, body, head, tail)
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
        swapV1 = chain(nBit, [put(nBit, (inBit,inBit+1)=>ConstGate.SWAP ) for irBit=rBit:-1:1 for inBit=irBit  :(vBit+irBit-1)])
        swap1V = chain(nBit, [put(nBit, (inBit+1,inBit)=>ConstGate.SWAP ) for irBit=1:   rBit for inBit=(vBit+irBit-1):-1:irBit])
        # println("s2\n")
        # dispatch!(cBlockHead, :random)
        cBlockHead = chain(nBit, DCbuilder(nBit, depth).block, swapV1) |> diff
        cBlocks = []
        cEBlocks = []        
        for i=1:nBlock
            # println("s3\n")
            cBlockHead = chain(nBit, DCbuilder(nBit, depth).block, swapV1) |> diff
            # println("s4\n")
            cBlock =  chain(nBit, vcat([swap1V], cBlockHead.blocks))
            push!(cBlocks, cBlock)
            # println("s4_1\n")
            push!(cEBlocks, put(nBitA, Tuple((nBitA-i*rBit-vBit+1):(nBitA-(i-1)*rBit),)=>cBlockHead))
            # println("s4_2\n")
        end
        cBlocks[1] = cEBlocks[1]
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
            push!(cBlocks, chain(nBit,put(nBit, (2,1)=>ConstGate.SWAP), put(nBit, 1=>H), control(nBit, 1, nBit=>Z)))
        end
        circuit = chain(nBit,cBlocks)  
        cExtend = chain(nBitA, vcat([repeat(nBitA,H,Tuple(nBitA:-1:1,))], [control(nBitA, i-1, i=>Z) for i=nBitA:-1:2])) 
    end
    
    (circuit, cExtend)
end