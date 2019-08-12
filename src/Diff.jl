#= 
Realizing functions of selecting differentiable blocks and doing differentiations.
Reference: https://github.com/QuantumBFS/QuAlgorithmZoo.jl.
=#
export Rotor, AbstractDiff, DiffBlock, QDiff, markDiff, getQdiff, getNdiff
import Yao: content, chcontent, mat, apply!
# Reminder of dependent packages.
## using StatsBase
## using MacroTools:@forward


# Basic type for differentiation.
abstract type AbstractDiff{GT, N, T} <: TagBlock{GT, N} end
Base.adjoint(df::AbstractDiff) = Daggered(df)
istraitkeeper(::AbstractDiff) = Val(true)
const Rotor{N, T} = Union{RotationGate{N, T}, PutBlock{N, <:Any, <:RotationGate{<:Any, T}}}
const CphaseGate{N, T} = ControlBlock{N,<:ShiftGate{T},<:Any}
const DiffBlock{N, T} = Union{Rotor{N, T}, CphaseGate{N, T}}


# Quantum differentiation block.
"""
    QDiff{GT, N, T} <: AbstractDiff{GT, N, T}
    QDiff(block) -> QDiff
Mark a block as quantum differentiable.
"""
mutable struct QDiff{GT, N, T} <: AbstractDiff{GT, N, T}
    block::GT
    grad::T
    QDiff(block::DiffBlock{N, T}) where {N, T} = new{typeof(block), N, T}(block, T(0))
end

content(cb::QDiff) = cb.block
chcontent(cb::QDiff, blk::DiffBlock) = QDiff(blk)
Base.adjoint(df::QDiff) = QDiff(content(df)')
@forward QDiff.block apply!
mat(::Type{T}, df::QDiff) where T = mat(T, df.block)

## Print the differentiation marks.
function YaoBlocks.print_annotation(io::IO, df::QDiff)
    printstyled(io, "[̂∂] "; bold=true, color=:yellow)
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
markDiff(block::Union{RotationGate, CphaseGate}) = QDiff(block)
## Exclude control blocks.
markDiff(block::ControlBlock) = block


"""
    getQdiff(psifunc::Function, diffblock::AbstractDiff, op::AbstractBlock) -> diffblock.grad::Float64
Quantum Operator differentiation.
\n `psifunc = ()-> reg::DefaultRegister |> c::ChainBlock`
\n `diffblock = collect_blocks(AbstractDiff, c)`
\n `op`: Witness Operator to measure reg.
"""
@inline function getQdiff(psifunc::Function, diffblock::AbstractDiff, op::AbstractBlock)
    r1, r2 = _perturb( ()->mean( expect(op, psifunc()) ) |> real, diffblock, π/2 )
    diffblock.grad = (r2 - r1)/2
end


"""
    getNdiff(psifunc::Function, parblock::AbstractBlock, op::AbstractBlock; δ::Real=0.01) -> diffblock.grad::Float64
Numerical Operator differentiation.
\n `psifunc = ()-> reg::DefaultRegister |> c::ChainBlock`
\n `diffblock`: Paramterized block(gate) in c.
\n `op`: Witness Operator to measure reg.
"""
@inline function getNdiff(psifunc::Function, parblock::AbstractBlock, op::AbstractBlock; δ::Real=0.01)
    r1, r2 = _perturb( ()->mean( expect(op, psifunc()) ) |> real, parblock, δ )
    parblock.grad = (r2 - r1) / (2δ)
end


@inline function _perturb(func, gate::AbstractDiff, δ::Real)
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end