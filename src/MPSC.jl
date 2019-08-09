#=
Main functions to build MPS-circuit.
=#
export MPSpar, MPSC, MPSDCpar


"""
    MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
Construct parameters that MPSC needs.
\nFields:
\n`nBitA::Int64`:  Number of lines(bits) in the Extended circuit. 
\n`vBit::Int64`:   Number of virtual bits in the MPS circuit.
\n`rBit::Int64`:   Number of reusable bits in the MPS circuit.
\n`nBit::Int64`:   Number of lines(bits) in the MPS circuit.
\n`nBlock::Int64`: Number of  MPS blocks in the MPS circuit. 
"""
struct MPSpar
    nBitA::Int64  # Number of lines(bits) in the Extended circuit. 
    vBit::Int64   # Number of virtual bits in the MPS circuit.
    rBit::Int64   # Number of reusable bits in the MPS circuit.
    nBit::Int64   # Number of lines(bits) in the MPS circuit.
    nBlock::Int64 # Number of  MPS blocks in the MPS circuit. 

    function MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
        if (nBitA - vBit -rBit) % rBit != 0 
            println("Error: nBlock is not integer!")
            return 0
        end
        nBlock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        new(nBitA, vBit, rBit, nBit, nBlock)
    end
end


"""
    MPSDCpar(circuit::ChainBlock)
Get the circuit parameters of a MPS-DC circuit or MPS-DC extended circuit.
\nFields:
\n`nBitA::Int64`:  Number of qubits(lines) in the MPS-DC extended circuit.
\n`vBit::Int64`:   Number of virtual qubits in the MPS-DC circuit.
\n`rBit::Int64`:   Number of reusable qubits in the MPS-DC circuit.
\n`nBit::Int64`:   Number of qubits(lines) in the MPS-DC circuit.
\n`nBlock::Int64`: Number of  MPS blocks in the MPS-DC circuit or MPS-DC extended circuit.
\n`depth::Int64`:  Depth(number of layers) of each DC block in the MPS-DC circuit.
"""
struct MPSDCpar
    nBitA::Int64  # Number of qubits(lines) in the MPS-DC extended circuit.
    vBit::Int64   # Number of virtual qubits in the MPS-DC circuit.
    rBit::Int64   # Number of reusable qubits in the MPS-DC circuit.
    nBit::Int64   # Number of qubits(lines) in the MPS-DC circuit.
    nBlock::Int64 # Number of  MPS blocks in the MPS-DC circuit or MPS-DC extended circuit.
    depth::Int64  # Depth(number of layers) of each DC block in the MPS-DC circuit.
    
    function MPSDCpar(circuit::ChainBlock)
        if typeof(circuit[1]) <: AbstractContainer
            # MPSC().cExtend[1]::AbstractContainer
            rBit = content(circuit[1])[2][1].locs[1]
            nBit = nqubits(content(circuit[1]))
            vBit = nBit - rBit
            depth = length(content(circuit[1])[1])
            nBitA = nqubits(circuit) 
        elseif typeof(circuit[1]) <: CompositeBlock
            # MPSC().circuit[1]::ChainBlock
            rBit = circuit[1][1][2][1].locs[1]
            nBit = nqubits(circuit)
            vBit = nBit - rBit
            depth = length(circuit[1][1][1])
            nBitA = vBit + rBit*length(circuit)
        else
            println("ERROR: Input circuit is not supported by the function!")
            return 1
        end
        nBlock = Int((nBitA - vBit) / rBit)
    new(nBitA, vBit, rBit, nBit, nBlock, depth)
    end
end


"""
    MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar=0)
Structure of related elements of MPS circuit.
\n`blockT` = ("DC", depth) stands for "Differentiable circuit".
\n`blockT` = "CS" stands for "cluster state".
\n`dBlocksPar` is the array of the parameters for Differentiable blocks in the circuit. 
\nFields:
\n`circuit::ChainBlock`:                MPS circuit.
\n`mpsBlocks::Array{CompositeBlock,1}`: Array of all the MPS blocks in the MPS circuit.
\n`cExtend::ChainBlock`:                The MPS circuit extended back to where it doesn't reuse any qubit.
\n`cEBlocks::Array{CompositeBlock,1}`:  Array of all the MPS blocks in the Extended circuit.
\n`dGates::Array{AbstractDiff,1}`:      Differentiable gates of the MPS circuit if applicable.
\n`nBit::Int64`:                        Number of lines(bits) of the MPS circuit. 
\n`nBlock::Int64`:                      Number of blocks in the MPS circuit.
"""
struct MPSC
    circuit::ChainBlock                # MPS circuit.
    mpsBlocks::Array{CompositeBlock,1} # Array of all the MPS blocks in the MPS circuit.
    cExtend::ChainBlock                # The MPS circuit extended back to where it doesn't reuse any qubit.
    cEBlocks::Array{CompositeBlock,1}  # Array of all the MPS blocks in the Extended circuit.
    dGates::Array{AbstractDiff,1}      # Differentiable gates of the MPS circuit if applicable.
    nBit::Int64                        # Number of lines(bits) of the MPS circuit. 
    nBlock::Int64                      # Number of blocks in the MPS circuit.

    function MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar::Array{Float64,1}=[0.0])
        par2nd = MPSpar(nBitA, vBit, rBit)
        nBlock = par2nd.nBlock
        nBit = par2nd.nBit
        MPS = MPSbuilder(nBitA, vBit, rBit, blockT)
        circuit = MPS.circuit
        cExtend = MPS.cExtend
        dispatch!(circuit, :random)
        dGates = collect_blocks(AbstractDiff, circuit)
        if dBlocksPar != [0.0]
            if length(dGates) == length(dBlocksPar)
                pars = [parameters(dGates[i])[1] for i=1:length(dGates)] 
                pars .= dBlocksPar
                dispatch!.(dGates, pars) 
            else
                println("The number of elements in dBlocksPar is not correct!!")
                return 1
            end
        end
        mpsBlocks = CompositeBlock[]
        for i = 1:length(circuit)-1
            push!(mpsBlocks, circuit[i][1])
        end 
        push!(mpsBlocks, circuit[end])
        new(circuit, mpsBlocks, cExtend, cExtend.blocks, dGates, nBit, nBlock)
    end
end