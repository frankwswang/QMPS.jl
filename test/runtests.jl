push!(LOAD_PATH, abspath("src"))
using MPSCircuit
using Test, Random, Statistics, Yao, Yao.ConstGate, YaoExtensions

@testset "CircuitBuilder.jl + MPSC.jl" begin
    #####CircuitBuilder.jl#####
    # DCbuilder(nBit::Int64, depth::Int64)
    seedNum = 1234
    dc = DCbuilder(4,2)

    head  = chain(4, chain(4, put(4,4=>Rx(0)), put(4,3=>Rx(0)), put(4,2=>Rx(0)), put(4,1=>Rx(0))),
                     chain(4, put(4,4=>Rz(0)), put(4,3=>Rz(0)), put(4,2=>Rz(0)), put(4,1=>Rz(0))),
                     chain(4, control(4,4,3=>X), control(4,3,2=>X), control(4,2,1=>X), control(4,1,4=>X)))
    block = chain(4, chain(4, put(4,4=>Rz(0)), put(4,3=>Rz(0)), put(4,2=>Rz(0)), put(4,1=>Rz(0))),
                     chain(4, put(4,4=>Rx(0)), put(4,3=>Rx(0)), put(4,2=>Rx(0)), put(4,1=>Rx(0))),
                     chain(4, put(4,4=>Rz(0)), put(4,3=>Rz(0)), put(4,2=>Rz(0)), put(4,1=>Rz(0))),
                     chain(4, control(4,4,3=>X), control(4,3,2=>X), control(4,2,1=>X), control(4,1,4=>X)))
    tail  = chain(4, chain(4, put(4,4=>Rz(0)), put(4,3=>Rz(0)), put(4,2=>Rz(0)), put(4,1=>Rz(0))),
                     chain(4, put(4,4=>Rx(0)), put(4,3=>Rx(0)), put(4,2=>Rx(0)), put(4,1=>Rx(0))))
    Cblock = chain(4, block, block)
    body = chain(4, head, block, tail)
    @test head == dc.head
    @test block == dc.block
    @test tail == dc.tail
    @test Cblock == dc.Cblock
    @test body == dc.body 

    # MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::String)
    regCSr = repeat(rand_state(2), 5000)
    regCSe = rand_state(4)
    mpsCS = MPSbuilder(4, 1, 1, "CS")
    CScr = mpsCS.circuit
    CSce = mpsCS.cExtend
    CScrt = chain(2, chain(2, repeat(2, H, (2,1)), control(2, 1, 2=>Z), Measure(2, locs=2, resetto=0)), 
                     chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z), Measure(2, locs=2, resetto=0)), 
                     chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z)))
    CScet = chain(4, repeat(4, H, (4,3,2,1)), control(4, 3, 4=>Z), control(4, 2, 3=>Z), control(4, 1, 2=>Z))
    opx(n) = chain(n, [put(n, i=>X) for i=1:n])
    opy(n) = chain(n, [put(n, i=>Y) for i=1:n])
    opz(n) = chain(n, [put(n, i=>Z) for i=1:n])
    function cr_test(reg::ArrayReg, cr::CompositeBlock, crt::CompositeBlock)
        n = nqubits(reg)
        Random.seed!(seedNum)
        reg1x = expect(opx(n), copy(reg) |> cr) |> mean |> real
        Random.seed!(seedNum)
        reg1y = expect(opy(n), copy(reg) |> cr) |> mean |> real
        Random.seed!(seedNum)
        reg1z = expect(opz(n), copy(reg) |> cr) |> mean |> real
        Random.seed!(seedNum)
        reg2x = expect(opx(n), copy(reg) |> crt) |> mean |> real
        Random.seed!(seedNum)
        reg2y = expect(opy(n), copy(reg) |> crt) |> mean |> real
        Random.seed!(seedNum)
        reg2z = expect(opz(n), copy(reg) |> crt) |> mean |> real
        @test reg1x ≈ reg2x
        @test reg1y ≈ reg2y
        @test reg1z ≈ reg2z
    end
    function ce_test(reg::ArrayReg, ce::CompositeBlock, cet::CompositeBlock)
        reg3 = copy(reg) |> ce
        reg4 = copy(reg) |> cet
        @test reg3 ≈ reg4
    end
    cr_test(regCSr, CScr, CScrt)
    ce_test(regCSe, CSce, CScet)

    # MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::Tuple{String, Int64})
    regDCr = repeat(rand_state(4), 5000)
    regDCe = rand_state(6)
    mpsDC = MPSbuilder(6, 2, 2, ("DC", 2))
    DCcr = mpsDC.circuit
    DCce = mpsDC.cExtend
    Cb1 = deepcopy(Cblock) |> markDiff
    Cb2 = deepcopy(Cblock) |> markDiff
    dispatch!(Cb1, :random)
    dispatch!(Cb2, :random)
    pars1 = getiparams.(content.(collect_blocks(QDiff, Cb1)))
    pars2 = getiparams.(content.(collect_blocks(QDiff, Cb2)))
    dispatch!(DCcr[1], pars1)
    dispatch!(DCcr[2], pars2)
    swap1 = chain(4, put(4, (2,3)=>SWAP), put(4, (3,4)=>SWAP), put(4, (1,2)=>SWAP), put(4, (2,3)=>SWAP))
    swap2 = chain(4, put(4, (3,2)=>SWAP), put(4, (2,1)=>SWAP), put(4, (4,3)=>SWAP), put(4, (3,2)=>SWAP))
    DCcrt = chain(4, Cb1, swap1, Measure(4, locs=(3,4), resetto=0), swap2, Cb2, swap1)
    DCcet = chain(6, subroutine(6, chain(4, Cb1, swap1), 3:6), subroutine(6, chain(4, Cb2, swap1), 1:4))
    cr_test(regDCr, DCcr, DCcrt)
    ce_test(regDCe, DCce, DCcet)

    #####MPSC.jl#####
    # MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
    nA = 13
    v = 4
    r = 3
    d = 2
    par = MPSpar(nA, v, r)
    @test par.nBitA == nA
    @test par.vBit == v
    @test par.rBit == r
    @test par.nBit == v+r
    @test par.nBlock == Int((nA - v) / r)

    # MPSDCpar(circuit::ChainBlock)
    mpsDC2 = MPSbuilder(nA, v, r, ("DC", d))
    par2 = MPSDCpar(mpsDC2.circuit)
    @test par2.nBitA == nA
    @test par2.vBit == v
    @test par2.rBit == r
    @test par2.nBit == v+r
    @test par2.nBlock == length(mpsDC2.circuit)
    @test par2.depth == d
    par3 = MPSDCpar(mpsDC2.cExtend)
    @test par3.nBitA == nA
    @test par3.vBit == v
    @test par3.rBit == r
    @test par3.nBit == v+r
    @test par3.nBlock == length(mpsDC2.cExtend)
    @test par3.depth == d

    # MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar=0)
    ## blockT = "CS"
    mpssCS = MPSC("CS", 4, 1, 1)
    cr_test(regCSr, mpssCS.circuit, CScrt)
    ce_test(regCSe, mpssCS.cExtend, CScet)
    @test mpssCS.mpsBlocks == [chain(2, repeat(2, H, (2,1)), control(2, 1, 2=>Z)), 
                               chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z)), 
                               chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z))]
    @test mpssCS.cEBlocks == CScet.blocks
    @test mpssCS.nBit == nqubits(CScrt)
    @test mpssCS.nBlock == length(mpssCS.circuit)                          
    ## blockT = ("DC", 2)
    mpssDC = MPSC(("DC", 2), 6, 2, 2)
    reg2 = rand_state(6)
    Random.seed!(seedNum)
    r1 = copy(reg2)|> mpssDC.cExtend
    Random.seed!(seedNum)
    r2 = copy(reg2)|> DCcet 
    @test (r1 ≈ r2) == false
    dG2 = collect_blocks(QDiff, mpssDC.circuit)
    dispatch!.(dG2, [pars1;pars2]) 
    cr_test(regDCr, mpssDC.circuit, DCcrt)
    ce_test(regDCe, mpssDC.cExtend, DCcet)
    @test mpssDC.mpsBlocks == [chain(4, Cb1, swap1), chain(4, swap2, Cb2, swap1)]
    @test mpssDC.cEBlocks == DCce.blocks
    @test mpssDC.dGates == collect_blocks(QDiff, DCce)
    @test mpssDC.nBit == nqubits(DCcr)
    @test mpssDC.nBlock == length(mpssDC.circuit)  
end

@testset "Diff.jl" begin
    # markDiff(block::AbstractBlock) QDiff(block)
    seedNum = 1234
    n = 3
    del = 10e-6 
    Random.seed!(seedNum)
    reg = rand_state(n, nbatch = 1000)
    c = chain(n, put(n, 1=>Rx(pi)), put(n, 1=>shift(0)), put(n, 2=>X), put(n, 1=>Ry(0)))
    @show c
    c = markDiff(c)
    @show c
    dG = collect_blocks(QDiff, c)
    @test dG == [QDiff(Rx(pi)), QDiff(Ry(0.0))]

    # getQdiff(psifunc, diffblock::QDiff, op::AbstractBlock)
    # getNdiff(overlapFunc::Function, dGate::QDiff; δ::Real=0.01)
    c2 = DCbuilder(n,3).body
    c2 = markDiff(c2)
    dispatch!(c2, :random)
    op = put(n, 2=>Z)
    dGates = collect_blocks(QDiff, c2)
    qg = getQdiff.(()->(copy(reg) |> c2), dGates, Ref(op))
    ng = getNdiff.(()->(copy(reg) |> c2), dGates, Ref(op), δ=del)
    tg = zeros(length(dGates))
    for i = 1:length(dGates)
        dispatch!(-, dGates[i], (del,))
        psi1 = expect(op, copy(reg) |> c2) |> mean |> real
        dispatch!(+, dGates[i], (2del,))
        psi2 = expect(op, copy(reg) |> c2) |> mean |> real
        dispatch!(-, dGates[i], (del,))
        tg[i] = (psi2 - psi1) / (2del)
    end
    @test tg ≈ qg
    @test ng ≈ qg
end