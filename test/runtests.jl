using QMPS
using Test, Yao
import Random: seed! 
import Statistics: mean 

@testset "CircuitBuilder.jl + MPSC.jl" begin
    #####Test Functions#####
    function cr_test(reg::ArrayReg, cr::CompositeBlock, crt::CompositeBlock)
        n = nqubits(reg)
        seed!(seedNum)
        ExpReg1x = expect(opx(n), copy(reg) |> cr) |> mean |> real
        seed!(seedNum)
        ExpReg1y = expect(opy(n), copy(reg) |> cr) |> mean |> real
        seed!(seedNum)
        ExpReg1z = expect(opz(n), copy(reg) |> cr) |> mean |> real
        seed!(seedNum)
        ExpReg2x = expect(opx(n), copy(reg) |> crt) |> mean |> real
        seed!(seedNum)
        ExpReg2y = expect(opy(n), copy(reg) |> crt) |> mean |> real
        seed!(seedNum)
        ExpReg2z = expect(opz(n), copy(reg) |> crt) |> mean |> real
        @test ExpReg1x ≈ ExpReg2x
        @test ExpReg1y ≈ ExpReg2y
        @test ExpReg1z ≈ ExpReg2z
    end

    function ce_test(reg::ArrayReg, ce::CompositeBlock, cet::CompositeBlock)
        reg3 = copy(reg) |> ce
        reg4 = copy(reg) |> cet
        @test reg3 ≈ reg4
    end

    function re_test(regr::ArrayReg, rege::ArrayReg, cr::CompositeBlock, ce::CompositeBlock, vBit::Int64)
        nr = nqubits(regr)
        ne = nqubits(rege)
        seed!(seedNum) 
        ExpRegrx = expect(opVx(nr, vBit), copy(regr) |> cr) |> mean |> real
        seed!(seedNum)
        ExpRegry = expect(opVy(nr, vBit), copy(regr) |> cr) |> mean |> real
        seed!(seedNum)
        ExpRegrz = expect(opVz(nr, vBit), copy(regr) |> cr) |> mean |> real
        seed!(seedNum)
        ExpRegex = expect(opVx(ne, vBit), copy(rege) |> ce) |> mean |> real
        seed!(seedNum)
        ExpRegey = expect(opVy(ne, vBit), copy(rege) |> ce) |> mean |> real
        seed!(seedNum)
        ExpRegez = expect(opVz(ne, vBit), copy(rege) |> ce) |> mean |> real
        @test isapprox(ExpRegrx, ExpRegex, atol = 0.5e-2)
        @test isapprox(ExpRegry, ExpRegey, atol = 0.5e-2)
        @test isapprox(ExpRegrz, ExpRegez, atol = 0.5e-2)
    end

    opx(n::Int64) = chain(n, [put(n, i=>X) for i=1:n])
    opy(n::Int64) = chain(n, [put(n, i=>Y) for i=1:n])
    opz(n::Int64) = chain(n, [put(n, i=>Z) for i=1:n])

    opVx(n::Int64, V::Int64) = chain(n, [put(n, i=>X) for i=1:V])
    opVy(n::Int64, V::Int64) = chain(n, [put(n, i=>Y) for i=1:V])
    opVz(n::Int64, V::Int64) = chain(n, [put(n, i=>Z) for i=1:V])


    #####CircuitBuilder.jl#####
    seedNum = 1234
    # DCbuilder(nBit::Int64, depth::Int64)
    ## Preparing section
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
    ## Testing section
    @test head == dc.head
    @test block == dc.block
    @test tail == dc.tail
    @test Cblock == dc.Cblock
    @test body == dc.body 

    # MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::String)
    ## Preparing section
    regCSr = rand_state(2, nbatch=5000)
    regCSe = rand_state(4)
    regCSr0 = zero_state(2, nbatch=10000)
    regCSe0 = zero_state(4, nbatch=10000)
    mpsCS = MPSbuilder(4, 1, 1, "CS") 
    CScr = mpsCS.circuit
    CSce = mpsCS.cExtend
    CScrt = chain(2, chain(2, repeat(2, H, (2,1)), control(2, 1, 2=>Z), Measure(2, locs=2, resetto=0)), 
                     chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z), Measure(2, locs=2, resetto=0)), 
                     chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z)))
    CScet = chain(4, repeat(4, H, (4,3,2,1)), control(4, 3, 4=>Z), control(4, 2, 3=>Z), control(4, 1, 2=>Z))
    ## Testing section
    print(IOBuffer(), mpsCS)
    re_test(regCSr0, regCSe0, CScr, CSce, 1)
    re_test(regCSr0, regCSe0, CScrt, CScet, 1)
    cr_test(regCSr, CScr, CScrt)
    ce_test(regCSe, CSce, CScet)
    
    # MPSbuilder(nBitA::Int64, vBit::Int64, rBit::Int64, blockT::Tuple{String, Int64})
    ## Preparing section
    regDCr = rand_state(4, nbatch=5000)
    regDCe = rand_state(6)
    regDCr0 = zero_state(4, nbatch=10000)
    regDCe0 = zero_state(6, nbatch=10000)
    mpsDC = MPSbuilder(6, 2, 2, ("DC", 2))
    DCcr = mpsDC.circuit
    DCce = mpsDC.cExtend
    Cb1 = deepcopy(Cblock) |> markDiff
    Cb2 = deepcopy(Cblock) |> markDiff
    dispatch!(Cb1, :random)
    dispatch!(Cb2, :random)
    pars1 = getiparams.(content.(collect_blocks(QMPS.QDiff, Cb1)))
    pars2 = getiparams.(content.(collect_blocks(QMPS.QDiff, Cb2)))
    dispatch!(DCcr[1], pars1)
    dispatch!(DCcr[2], pars2)
    swap1 = chain(4, put(4, (2,3)=>SWAP), put(4, (3,4)=>SWAP), put(4, (1,2)=>SWAP), put(4, (2,3)=>SWAP))
    swap2 = chain(4, put(4, (3,2)=>SWAP), put(4, (2,1)=>SWAP), put(4, (4,3)=>SWAP), put(4, (3,2)=>SWAP))
    DCcrt = chain(4, Cb1, swap1, Measure(4, locs=(3,4), resetto=0), swap2, Cb2, swap1)
    DCcet = chain(6, subroutine(6, chain(4, Cb1, swap1), 3:6), subroutine(6, chain(4, Cb2, swap1), 1:4))
    ## Testing section
    print(IOBuffer(), mpsDC)
    re_test(regDCr0, regDCe0, DCcr, DCce, 2)
    re_test(regDCr0, regDCe0, DCcrt, DCcet, 2)
    cr_test(regDCr, DCcr, DCcrt)
    ce_test(regDCe, DCce, DCcet)


    #####MPSC.jl#####
    # MPSpar(nBitA::Int64, vBit::Int64, rBit::Int64)
    ## Preparing section
    nA = 13
    v = 4
    r = 3
    d = 2
    par = MPSpar(nA, v, r)
    ## Testing section
    @test par.nBitA == nA
    @test par.vBit == v
    @test par.rBit == r
    @test par.nBit == v+r
    @test par.nBlock == Int((nA - v) / r)

    # MPSDCpar(circuit::ChainBlock)
    ## Preparing section
    mpsDC2 = MPSbuilder(nA, v, r, ("DC", d))
    par2 = MPSDCpar(mpsDC2.circuit)
    par3 = MPSDCpar(mpsDC2.cExtend)
    ## Testing section
    @test par2.nBitA == nA
    @test par2.vBit == v
    @test par2.rBit == r
    @test par2.nBit == v+r
    @test par2.nBlock == length(mpsDC2.circuit)
    @test par2.depth == d
    @test par3.nBitA == nA
    @test par3.vBit == v
    @test par3.rBit == r
    @test par3.nBit == v+r
    @test par3.nBlock == length(mpsDC2.cExtend)
    @test par3.depth == d

    # MPSC(blockT, nBitA::Int64, vBit::Int64, rBit::Int64=1; dBlocksPar=0)
    ## blockT = "CS"
    ### Preparing section
    mpssCS = MPSC("CS", 4, 1, 1)
    ### Testing section
    re_test(regCSr0, regCSe0, mpssCS.circuit, mpssCS.cExtend, 1)
    cr_test(regCSr, mpssCS.circuit, CScrt)
    ce_test(regCSe, mpssCS.cExtend, CScet)
    @test mpssCS.mpsBlocks == [chain(2, repeat(2, H, (2,1)), control(2, 1, 2=>Z)), 
                               chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z)), 
                               chain(2, put(2, (2,1)=>SWAP), put(2, 1=>H), control(2, 1, 2=>Z))]
    @test mpssCS.cEBlocks == CScet.blocks
    @test mpssCS.nBit == nqubits(CScrt)
    @test mpssCS.nBlock == length(mpssCS.circuit)
    ## blockT = ("DC", 2)
    ### Preparing section
    mpssDC = MPSC(("DC", 2), 6, 2, 2)
    reg2DCe = rand_state(6)
    seed!(seedNum)
    reg2DCe_1 = copy(reg2DCe)|> mpssDC.cExtend
    seed!(seedNum)
    reg2DCe_2 = copy(reg2DCe)|> DCcet 
    seed!(seedNum)
    regDCr0_1 = copy(regDCr0)|> mpssDC.circuit
    seed!(seedNum)
    regDCr0_2 = copy(regDCr0)|> DCcrt
    dG2 = collect_blocks(QMPS.QDiff, mpssDC.circuit)
    dispatch!.(dG2, [pars1;pars2])
    ### Testing section
    @test (reg2DCe_1 ≈ reg2DCe_2) == false
    @test (regDCr0_1 ≈ regDCr0_2) == false
    cr_test(regDCr0, mpssDC.circuit, DCcrt)
    ce_test(reg2DCe, mpssDC.cExtend, DCcet)
    re_test(regDCr0, regDCe0, mpssDC.circuit, mpssDC.cExtend, 2)
    @test mpssDC.mpsBlocks == [chain(4, Cb1, swap1), chain(4, swap2, Cb2, swap1)]
    @test mpssDC.cEBlocks == DCce.blocks
    @test mpssDC.dGates == collect_blocks(QMPS.QDiff, DCce)
    @test mpssDC.nBit == nqubits(DCcr)
    @test mpssDC.nBlock == length(mpssDC.circuit)
    @test parameters(MPSC(("DC", 1), 3, 1, 1, dBlocksPar=[1.0:12.0;]).circuit) == [1.0:12.0;]  


    #####Error Exception Tests#####
    @test_throws ErrorException MPSbuilder(6, 3, 2, ("DC", 2))
    @test_throws ErrorException MPSpar(5, 1, 3)
    @test_throws ErrorException mpssDC2 = MPSC(("DC", 1), 5, 2, 2)

    @test_throws ErrorException MPSbuilder(4, 2, 1, "CS") 
    @test_throws ErrorException MPSbuilder(5, 1, 2, "CS")

    @test_throws ErrorException MPSDCpar(chain(chain(4, Cb1, put(4, 1=>X))))
    @test_throws ErrorException mpssDC2 = MPSC("CS", 3, 1, 1, dBlocksPar=[1:0.5:13;])
    @test_throws ErrorException mpssDC2 = MPSC(("DC", 1), 3, 1, 1, dBlocksPar=[1:0.5:13;])
end

@testset "Diff.jl" begin
    # Preparing section
    seedNum = 1234
    n = 3
    
    # Testing basic Yao-compatible functions for QDiff{GT, N} and markDiff(block::AbstractBlock).
    c = chain(n, put(n, 1=>Rx(pi)), put(n, 1=>shift(0)), put(n, 2=>X), control(n, 2, 1=>Ry(0)), control(n, 2, 1=>shift(pi/2)))
    c = markDiff(c)
    DB_Rx = collect_blocks(QMPS.QDiff, c)[1]
    DB_CS = collect_blocks(QMPS.QDiff, c)[end]
    applyTestReg1 = rand_state(1)
    applyTestReg2 = rand_state(n)
    @test (copy(applyTestReg1) |> DB_Rx) == (apply!(copy(applyTestReg1), DB_Rx))
    @test (copy(applyTestReg1) |> DB_Rx.block) == (apply!(copy(applyTestReg1), DB_Rx.block))
    @test (copy(applyTestReg1) |> DB_Rx) == (copy(applyTestReg1) |> DB_Rx.block)
    @test (copy(applyTestReg2) |> DB_CS) == (apply!(copy(applyTestReg2), DB_CS))
    @test (copy(applyTestReg2) |> DB_CS.block) == (apply!(copy(applyTestReg2), DB_CS.block))
    @test (copy(applyTestReg2) |> DB_CS) == (copy(applyTestReg2) |> DB_CS.block)
    dG = collect_blocks(QMPS.QDiff, c)
    @test dG == [ QMPS.QDiff(Rx(pi)), QMPS.QDiff(control(n, 2, 1=>shift(pi/2))) ]
    @test dG[1].mat == mat(dG[1])
    @test (dG[1])' == adjoint(dG[1]) == QMPS.QDiff( adjoint(dG[1].block) )
    
    # Testing functions getQdiff!() and getNdiff().
    del = 10e-6 
    seed!(seedNum)
    reg = rand_state(n, nbatch = 1000)
    c2 = DCbuilder(n,3).body
    c2 = markDiff(c2)
    dispatch!(c2, :random)
    op = put(n, 2=>Z)
    dGates = collect_blocks(QMPS.QDiff, c2)
    qg = getQdiff!.(()->(copy(reg) |> c2), dGates, Ref(op))
    ng = getNdiff.(()->(copy(reg) |> c2), dGates, Ref(op), δ=del)
    qg_m2 = getQdiff!(()->(copy(reg) |> c2), dGates, op)
    ng_m2 = getNdiff(()->(copy(reg) |> c2), dGates, op, δ=del)
    tg = zeros(length(dGates))
    for i = 1:length(dGates)
        dispatch!(-, dGates[i], (del,))
        psi1 = expect(op, copy(reg) |> c2) |> mean |> real
        dispatch!(+, dGates[i], (2del,))
        psi2 = expect(op, copy(reg) |> c2) |> mean |> real
        dispatch!(-, dGates[i], (del,))
        tg[i] = (psi2 - psi1) / (2del)
    end
    @test qg == qg_m2
    @test ng == ng_m2
    @test tg ≈ qg
    @test ng ≈ qg
end