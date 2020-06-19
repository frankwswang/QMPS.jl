export markDiff, getQdiff, getNdiff
import StatsBase: mean


#### Structures #####
# Define supertypes of differentiable blocks.
const CphaseGate{N, T} = ControlBlock{N, <:ShiftGate{T}, <:Any}
const DiffBlock{N, T} = Union{RotationGate{N, T, <:Any}, CphaseGate{N, T}}


# Define the struct of a differentiable gate `QDiff{GT, N}`.
"""
    QDiff{GT, N} <: TagBlock{GT, N}
Differentiable gate. If you want to mark or return a differentiable block, use `markDiff` or `markDiff!` instead.
\nFields:
\n`block::GT`: Sub-block(s) of the differentiable block.  
\n`grad::Float64`: Gradient(s) of the differentiable block.
"""
mutable struct QDiff{GT, N} <: TagBlock{GT, N}
    block::GT
    grad::Float64
    mat
    QDiff(blk::DiffBlock{N}) where {N} = new{typeof(blk), N}(blk, 0, mat(blk))
end


#### Functions ##### 
# Make Marks of differentiable gates when printing circuit layout.
function YaoBlocks.print_annotation(io::IO, df::QDiff)
    printstyled(io, "[̂∂(MPSC)] "; bold=true, color=:yellow)
end


# Define basic Yao-compatible functions for QDiff{GT, N}.
Yao.content(cb::QDiff) = cb.block
Yao.apply!(reg::AbstractRegister, db::QDiff) = apply!(reg, Yao.content(db))
Yao.mat(::Type{T}, df::QDiff) where T = Yao.mat(T, df.block)
Base.adjoint(df::QDiff) = QDiff(Yao.content(df)')


# Convert gates that are differentiable into type `QDiff{GT, N}`.
"""
    markDiff(block::AbstractBlock) -> block::AbstractBlock
    Return the differentiable gate(s) `QDiff{GT, N}` from a block or a block tree such as `ChainBlock`.
"""
## For Block trees.
function markDiff(blk::AbstractBlock)
    blks = subblocks(blk)
    isempty(blks) ? blk : chsubblocks(blk, markDiff.(blks))
end
## For differentiable blocks.
markDiff(block::DiffBlock) = QDiff(block)
## Exclude control blocks.
markDiff(block::ControlBlock) = block


"""
    getQdiff(psifunc::Function, diffblock::QMPS.QDiff, op::AbstractBlock) -> diffblock.grad::Float64
Quantum Operator differentiation.
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `diffblock = collect_blocks(QMPS.QDiff, c)`
\n `op`: Witness Operator to measure reg.
"""
@inline function getQdiff(psifunc::Function, diffblock::QDiff, op::AbstractBlock)
    r1, r2 = _perturb( ()->mean( Yao.expect(op, psifunc()) ) |> real, diffblock, π/2 )
    diffblock.grad = (r2 - r1)/2
end


"""
    getNdiff(psifunc::Function, parblock::QMPS.QDiff, op::AbstractBlock; δ::Real=0.01) -> parblock.grad::Float64
Numerical Operator differentiation.
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `diffblock`: Paramterized block(gate) in c.
\n `op`: Witness Operator to measure reg.
"""
@inline function getNdiff(psifunc::Function, parblock::QDiff, op::AbstractBlock; δ::Real=0.01)
    r1, r2 = _perturb( ()->mean( Yao.expect(op, psifunc()) ) |> real, parblock, δ )
    parblock.grad = (r2 - r1) / (2δ)
end


@inline function _perturb(func, gate::QDiff, δ::Real)
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end