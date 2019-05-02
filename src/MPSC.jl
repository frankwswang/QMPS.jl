"""
Main functions to build MPS-circuit.
"""

export setMPSpar, MPSC

# Structure of parameters that MPSC needs.
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



# Structure of elements MPSCircuit.jl provides.
struct MPSC
    circuit  # MPS circuit.
    cBlocks  # Array of all the MPS blocks in the MPS circuit.
    cExtend  # The MPS circuit extended back to where it doesn't reuse any qubit.
    cEBlocks # Array of all the MPS blocks in the Extended circuit.
    dGates   # Differentiable gates of the MPS circuit if applicable.
    nBit     # Number of lines(bits) of the MPS circuit. 
    nBlock   # Number of blocks in the MPS ciruict.

    function MPSC(MPSblock, nBitA::Int64, vBit::Int64, rBit::Int64=1;dBlocksPar=0)
        par2nd = setMPSpar(nBitA, vBit, rBit)
        nBlock = par2nd.nBlock
        nBit = par2nd.nBit
        #println("M1\n")
        MPS = MPSbuilder(nBitA, vBit, rBit, MPSblock)
        #println("M2\n")
        circuit = MPS[1]
        #println("M3\n")
        cExtend = MPS[2]
        dGates = collect_blocks(AbstractDiff, circuit)
        if typeof(dBlocksPar) <: Array
            if length(dGates) == length(dBlocksPar)
                pars = [parameters(dGates[i])[1] for i=1:length(dGates)] 
                pars .= dBlocksPar
                dispatch!.(dGates, pars) 
            else
                println("The number of elements in dBlocksPar is not correct!!")
                return 1
            end
        elseif dBlocksPar != 0
            println("Invalid parameters input!! parameters should be collected in Arrays.")
                return 1
        else
        end 
        new(circuit, circuit.blocks, cExtend, cExtend.blocks, dGates, nBit, nBlock)
    end

end