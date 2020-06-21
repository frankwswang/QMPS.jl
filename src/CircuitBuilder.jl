#=
Circuit builders for MPSCircuit.jl.
=#
export DCbuilder, MPSbuilder


"""
    DCbuilder(nBit::Int64, depth::Int64)
Structure of elements may needed for a Quantum differentiable circuit.
\n
\nFields:
\n`block::ChainBlock`:    The block for 1 depth.
\n`circuit::ChainBlock`:  The standard structure of a differentiable circuit without the modifications on the head and the tail.
\n`fullbody::ChainBlock`: The modified version of a differentiable circuit with structures of a head and a tail (same number of blocks).
\n`head::ChainBlock`:     The "head" of the differentiable circuit(if directly input bits).
\n`tail::ChainBlock`:     The "tail" of the differentiable circuit(2 layers of rotaional gates). 
"""
struct DCbuilder
    block::ChainBlock    # The block for 1 depth.
    circuit::ChainBlock  # The standard structure of a differentiable circuit without the modifications on the head and the tail..
    fullbody::ChainBlock # The modified version of a differentiable circuit with structures of a head and a tail (same number of blocks). 
    head::ChainBlock     # The "head" of the differentiable circuit(if directly input bits).
    tail::ChainBlock     # The "tail" of the differentiable circuit(2 layers of rotaional gates).

    function DCbuilder(nBit::Int64, depth::Int64) #DON'T USE repeat to avoid POINTER side-effect!
        block = chain(nBit, vcat([chain(nBit,[put(nBit, i=>Rz(0)) for i=nBit:-1:1])], 
                                 [chain(nBit,[put(nBit, i=>Rx(0)) for i=nBit:-1:1])], 
                                 [chain(nBit,[put(nBit, i=>Rz(0)) for i=nBit:-1:1])], 
                                 [chain(nBit, vcat([control(nBit, i, (i-1)=>X) for i=nBit:-1:2], 
                                                   [control(nBit, 1, nBit=>X)]))] ))
        tail = chain(nBit, vcat([chain(nBit,[put(nBit, i=>Rz(0)) for i=nBit:-1:1])], 
                                [chain(nBit,[put(nBit, i=>Rx(0)) for i=nBit:-1:1])]))
        head = deepcopy(block[2:end])
        blocks = ChainBlock[]
        for i=1:depth push!(blocks, deepcopy(block)) end
        Cblock = chain(nBit, blocks)
        body = chain(nBit, vcat([head], blocks[2:end], [tail])) 
        new(block, Cblock, body, head, tail)
    end
end


"""
    MPSc{circuit::ChainBlock, cExtend::ChainBlock}
\n
Fields:
\n`circuit::ChainBlock`: QMPS circuit.
\n`cExtend::ChainBlock`: QMPS-extended circuit.
"""
struct MPSc
    circuit::ChainBlock # MPS-qubit-reusable circuit. 
    cExtend::ChainBlock # MPS-extended circuit.
end


"""
    MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::Union{String, Tuple{String, Int64}}) 
    -> 
    MPSc{circuit::ChainBlock, cExtend::ChainBlock}
Function for creating different types of MPS circuits.
\n`blockT::Union{String, Tuple{String, Int64}}`
\n1) `blockT = ("DC", depth)`
\n   "DC" stands for "differentiable circuit".
\n2) `blockT = "CS"`         
\n   "CS" stands for "cluster state".
\n`dBlocksPar::Array{Float64,1}` 
\n   The array of the parameters (for differentiable blocks) in the circuit.
"""
function MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::Union{String, Tuple{String, Int64}}) 
    par2nd = MPSpar(nBitA, vBit, rBit)
    nBlock = par2nd.nBlock
    nBit = par2nd.nBit
    if typeof(blockT) == Tuple{String, Int64}
        if blockT[1] == "DC" # differentiable circuit.
            depth = blockT[2]
            swapV1 = chain(nBit, [put(nBit, (inBit,inBit+1)=>SWAP) for irBit=rBit:-1:1 for inBit=irBit  :(vBit+irBit-1)])
            swap1V = chain(nBit, [put(nBit, (inBit+1,inBit)=>SWAP) for irBit=1:   rBit for inBit=(vBit+irBit-1):-1:irBit])
            MeasureBlock = Measure(nBit, locs=(nBit-vBit+1):nBit, resetto=0)
            cBlockHead = chain(nBit, chain(nBit, DCbuilder(nBit, depth).circuit, swapV1) |> markDiff, MeasureBlock)
            cBlocks = [cBlockHead]
            cEBlocks = [subroutine( nBitA, cBlockHead[1], (nBitA-rBit-vBit+1):nBitA )]      
            for i=2:nBlock-1
                cBlockHead = chain(nBit, DCbuilder(nBit, depth).circuit, swapV1) |> markDiff 
                cBlock =  chain(nBit, chain(nBit, vcat([swap1V], cBlockHead.blocks)), MeasureBlock)
                push!(cBlocks, cBlock)
                # push!(cEBlocks, put(nBitA, Tuple((nBitA-i*rBit-vBit+1):(nBitA-(i-1)*rBit),)=>cBlockHead)) #Use `subroutine` instead of `put` for CuYao compatibility and efficiency.
                push!(cEBlocks, subroutine( nBitA, cBlockHead, (nBitA-i*rBit-vBit+1):(nBitA-(i-1)*rBit) ))
            end
            cBlockHead = chain(nBit, DCbuilder(nBit, depth).circuit, swapV1) |> markDiff 
            cBlock =  chain(nBit, vcat([swap1V], cBlockHead.blocks))
            push!(cBlocks, cBlock)
            push!(cEBlocks, subroutine( nBitA, cBlockHead, (nBitA-nBlock*rBit-vBit+1):(nBitA-(nBlock-1)*rBit) ))
            circuit = chain(nBit, cBlocks)
            cExtend = chain(nBitA, cEBlocks)
        end
    end
    if typeof(blockT) == String
        if blockT == "CS" # MPS blocks for 1D cluster state.
            if vBit !=1 || rBit !=1
                error("The nBit=vBit+rBit of cluster state MPS blocks should be 2.\n")
            end
            cBlocks = ChainBlock[]
            MeasureBlock = Measure(nBit, locs=(nBit-vBit+1):nBit, resetto=0)
            push!(cBlocks, chain(nBit, chain(nBit, repeat(nBit,H,(2,1)), control(nBit, 1, nBit=>Z)), MeasureBlock))
            for i =2:nBlock-1
                push!(cBlocks, chain(nBit, chain(nBit,put(nBit, (2,1)=>SWAP), put(nBit, 1=>H), control(nBit, 1, nBit=>Z)), MeasureBlock))
            end
            push!(cBlocks, chain(nBit,put(nBit, (2,1)=>SWAP), put(nBit, 1=>H), control(nBit, 1, nBit=>Z)))
            circuit = chain(nBit,cBlocks)  
            cExtend = chain(nBitA, vcat([repeat(nBitA,H,Tuple(nBitA:-1:1,))], [control(nBitA, i-1, i=>Z) for i=nBitA:-1:2])) 
        end
    end

    MPSc(circuit, cExtend)
end