export QDiff, markDiff, getQdiff, getNdiff
import Yao: content, chcontent, mat, apply!
# Reminder of dependent packages.
## using StatsBase


# Basic type for differentiation.
const CphaseGate{N, T} = ControlBlock{N,<:ShiftGate{T},<:Any}
const DiffBlock{N, T} = Union{RotationGate{N, T, <:Any}, CphaseGate{N, T}}

# Quantum differentiation block.
"""
    QDiff{GT, N, T} <: TagBlock{GT, N}
    QDiff(block) -> QDiff
Mark a block as quantum differentiable.
"""
mutable struct QDiff{GT, N} <: TagBlock{GT, N}
    block::GT
    grad::Float64
    QDiff(block::DiffBlock{N}) where {N} = new{typeof(block), N}(block, 0)
end

content(cb::QDiff) = cb.block
chcontent(cb::QDiff, blk::AbstractBlock) = QDiff(blk)
YaoBlocks.PropertyTrait(::QDiff) = Yao.PreserveAll()

apply!(reg::AbstractRegister, db::QDiff) = apply!(reg, content(db))
mat(::Type{T}, df::QDiff) where T = mat(T, df.block)
Base.adjoint(df::QDiff) = QDiff(content(df)')

## Print the differentiation marks.
function YaoBlocks.print_annotation(io::IO, df::QDiff)
    printstyled(io, "[̂∂(MPSC)] "; bold=true, color=:yellow)
end


#### Functions ##### 
"""
    markDiff(block::AbstractBlock) -> block::AbstractBlock
Mark differentiable items in a block tree(e.g.: ChainBlock) as differentiable.
"""
function markDiff(blk::AbstractBlock)
    blks = subblocks(blk)
    isempty(blks) ? blk : chsubblocks(blk, markDiff.(blks))
end
## For differentiable blocks.
markDiff(block::DiffBlock) = QDiff(block)
## Exclude control blocks.
markDiff(block::ControlBlock) = block


"""
    getQdiff(psifunc::Function, diffblock::QDiff, op::AbstractBlock) -> diffblock.grad::Float64
Quantum Operator differentiation.
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `diffblock = collect_blocks(QDiff, c)`
\n `op`: Witness Operator to measure reg.
"""
@inline function getQdiff(psifunc::Function, diffblock::QDiff, op::AbstractBlock)
    r1, r2 = _perturb( ()->mean( expect(op, psifunc()) ) |> real, diffblock, π/2 )
    diffblock.grad = (r2 - r1)/2
end


"""
    getNdiff(psifunc::Function, parblock::AbstractBlock, op::AbstractBlock; δ::Real=0.01) -> diffblock.grad::Float64
Numerical Operator differentiation.
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `diffblock`: Paramterized block(gate) in c.
\n `op`: Witness Operator to measure reg.
"""
@inline function getNdiff(psifunc::Function, parblock::QDiff, op::AbstractBlock; δ::Real=0.01)
    r1, r2 = _perturb( ()->mean( expect(op, psifunc()) ) |> real, parblock, δ )
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