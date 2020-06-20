#=
Main functions to build MPS-circuit.
=#
export MPSpar, MPSC, MPSDCpar


"""
    MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
Construct parameters that MPSC needs.
\n
\nFields:
\n`nBitA::Int64`:  Number of lines(bits) in the QMPS Extended circuit.
\n`vBit::Int64`:   Number of virtual bits in the QMPS circuit.
\n`rBit::Int64`:   Number of reusable bits in the QMPS circuit.
\n`nBit::Int64`:   Number of lines(bits) in the QMPS circuit.
\n`nBlock::Int64`: Number of  MPS blocks in the QMPS circuit.
"""
struct MPSpar
    nBitA::Int64  # Number of lines(bits) in the QMPS Extended circuit.
    vBit::Int64   # Number of virtual bits in the QMPS circuit.
    rBit::Int64   # Number of reusable bits in the QMPS circuit.
    nBit::Int64   # Number of lines(bits) in the QMPS circuit.
    nBlock::Int64 # Number of  MPS blocks in the QMPS circuit.

    function MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
        if (nBitA - vBit -rBit) % rBit != 0
            error("Error: nBlock is not an integer.")
        end
        nBlock = Int((nBitA - vBit) / rBit)
        nBit = rBit + vBit
        new(nBitA, vBit, rBit, nBit, nBlock)
    end
end


"""
    MPSDCpar(circuit::ChainBlock)
Get the circuit parameters of a differentiable QMPS circuit (QMPS-DC) or of a QMPS-DC extended circuit.
\n
\nFields:
\n`nBitA::Int64`:  Number of qubits(lines) in the QMPS-DC extended circuit.
\n`vBit::Int64`:   Number of virtual qubits in the circuit.
\n`rBit::Int64`:   Number of reusable qubits in the circuit.
\n`nBit::Int64`:   Number of qubits(lines) in the circuit.
\n`nBlock::Int64`: Number of  MPS blocks in the circuit.
\n`depth::Int64`:  Depth(number of layers) of each DC block in the circuit.
"""
struct MPSDCpar
    nBitA::Int64  # Number of qubits(lines) in the MPS-DC extended circuit.
    vBit::Int64   # Number of virtual qubits in the MPS-DC circuit.
    rBit::Int64   # Number of reusable qubits in the MPS-DC circuit.
    nBit::Int64   # Number of qubits(lines) in the MPS-DC circuit.
    nBlock::Int64 # Number of  MPS blocks in the MPS-DC circuit or MPS-DC extended circuit.
    depth::Int64  # Depth(number of layers) of each DC block in the MPS-DC circuit.

    function MPSDCpar(circuit::ChainBlock)
        rBit = 0
        nBit = 0
        vBit = 0
        depth = 0
        nBitA = 0
        nBlock = 0
        try
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
            end
            nBlock = Int((nBitA - vBit) / rBit)    
        catch err
            if isa(err, DomainError) || (0 in (nBitA, vBit, rBit, nBit, nBlock, depth))
                error("ERROR: Input circuit is not supported by the function!")
            end
        end
        new(nBitA, vBit, rBit, nBit, nBlock, depth)
    end
end


"""
    MPSC(blockT::Union{String, Tuple{String, Int64}}, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar::Array{Float64,1}=[0.0])
Structure of elements related to a QMPS circuit.
\n`blockT::Union{String, Tuple{String, Int64}}`
\n1) `blockT = ("DC", depth)`
\n   "DC" stands for "differentiable circuit".
\n2) `blockT = "CS"`         
\n   "CS" stands for "cluster state".
\n`dBlocksPar::Array{Float64,1}` 
\n   The array of the parameters (for differentiable blocks) in the circuit.
\n
\nFields:
\n`circuit::ChainBlock`:                QMPS circuit.
\n`mpsBlocks::Array{CompositeBlock,1}`: Array of all the MPS blocks in the QMPS circuit.
\n`cExtend::ChainBlock`:                The circuit QMPS circuit is extended back to that doesn't reuse any qubit.
\n`cEBlocks::Array{CompositeBlock,1}`:  Array of all the MPS blocks in the Extended circuit.
\n`dGates::Array{QMPS.QDiff,1}`:        Differentiable gates of the QMPS circuit if applicable.
\n`nBit::Int64`:                        Number of lines (bits) of the QMPS circuit.
\n`nBlock::Int64`:                      Number of blocks in the QMPS circuit.
"""
struct MPSC
    circuit::ChainBlock                # QMPS circuit.
    mpsBlocks::Array{CompositeBlock,1} # Array of all the MPS blocks in the QMPS circuit.
    cExtend::ChainBlock                # The circuit QMPS circuit is extended back to that doesn't reuse any qubit.
    cEBlocks::Array{CompositeBlock,1}  # Array of all the MPS blocks in the Extended circuit.
    dGates::Array{QMPS.QDiff,1}        # Differentiable gates of the QMPS circuit if applicable.
    nBit::Int64                        # Number of lines (bits) of the QMPS circuit.
    nBlock::Int64                      # Number of blocks in the QMPS circuit.

    function MPSC(blockT::Union{String, Tuple{String, Int64}}, nBitA::Int64, vBit::Int64, rBit::Int64; dBlocksPar::Array{Float64,1}=[0.0])
        par2nd = MPSpar(nBitA, vBit, rBit)
        nBlock = par2nd.nBlock
        nBit = par2nd.nBit
        MPS = MPSbuilder(nBitA, vBit, rBit, blockT)
        circuit = MPS.circuit
        cExtend = MPS.cExtend
        dispatch!(circuit, :random)
        dGates = collect_blocks(QMPS.QDiff, circuit)
        if dBlocksPar != [0.0]
            if length(dGates) == length(dBlocksPar)
                pars = [parameters(dGates[i])[1] for i=1:length(dGates)]
                pars .= dBlocksPar
                dispatch!.(dGates, pars)
            elseif blockT == "CS"
                error("Cluster states don't support optional argument `dBlocksPar`.")
            else
                error("The number of elements in `dBlocksPar` is not correct.")
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