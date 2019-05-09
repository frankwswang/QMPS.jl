#=
Main functions to build MPS-circuit.
=#
export setMPSpar, MPSC


"""
    setMPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
Construct parameters that MPSC needs.
\nFields:
\n`nBitA::Int64`:  Number of lines(bits) in the Extended circuit. 
\n`vBit::Int64`:   Number of virtual bits in the MPS circuit.
\n`rBit::Int64`:   Number of reusable bits in the MPS circuit.
\n`nBit::Int64`:   Number of lines(bits) in the MPS circuit.
\n`nBlock::Int64`: Number of  MPS blocks in the MPS circuit. 
"""
struct setMPSpar
    nBitA::Int64  # Number of lines(bits) in the Extended circuit. 
    vBit::Int64   # Number of virtual bits in the MPS circuit.
    rBit::Int64   # Number of reusable bits in the MPS circuit.
    nBit::Int64   # Number of lines(bits) in the MPS circuit.
    nBlock::Int64 # Number of  MPS blocks in the MPS circuit. 

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


"""
    MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar=0)
Structure of related elements of MPS circuit.
\n`blockT` = ("DC", depth) stands for "Differentiable circuit".
\n`blockT` = "CS" stands for "cluster state".
\n `dBlocksPar` is the array of the parameters for Differentiable blocks in the circuit. 
\nFields:
\n`circuit::ChainBlock`:               MPS circuit.
\n`cBlocks::Array{CompositeBlock,1}`:  Array of all the MPS blocks in the MPS circuit.
\n`cExtend::ChainBlock`:               The MPS circuit extended back to where it doesn't reuse any qubit.
\n`cEBlocks::Array{CompositeBlock,1}`: Array of all the MPS blocks in the Extended circuit.
\n`dGates::Array{AbstractDiff,1}`:     Differentiable gates of the MPS circuit if applicable.
\n`nBit::Int64`:                       Number of lines(bits) of the MPS circuit. 
\n`nBlock::Int64`:                     Number of blocks in the MPS ciruict.

"""
struct MPSC
    circuit::ChainBlock               # MPS circuit.
    cBlocks::Array{CompositeBlock,1}  # Array of all the MPS blocks in the MPS circuit.
    cExtend::ChainBlock               # The MPS circuit extended back to where it doesn't reuse any qubit.
    cEBlocks::Array{CompositeBlock,1} # Array of all the MPS blocks in the Extended circuit.
    dGates::Array{AbstractDiff,1}     # Differentiable gates of the MPS circuit if applicable.
    nBit::Int64                       # Number of lines(bits) of the MPS circuit. 
    nBlock::Int64                     # Number of blocks in the MPS ciruict.

    function MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar::Array{Float64,1}=[0.0])
        par2nd = setMPSpar(nBitA, vBit, rBit)
        nBlock = par2nd.nBlock
        nBit = par2nd.nBit
        MPS = MPSbuilder(nBitA, vBit, rBit, blockT)
        circuit = MPS[1]
        cExtend = MPS[2]
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
        new(circuit, circuit.blocks, cExtend, cExtend.blocks, dGates, nBit, nBlock)
    end
end