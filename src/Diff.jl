"""
Realizing functions of selecting differentiable blocks and doing differentiations.
"""

export Rotor, AbstractDiff, DiffBlock, QDiff, diff, opdiff
import Yao: expect, content, chcontent, mat, apply!
# Reminder of dependent packages.
## using StatsBase
## using MacroTools:@forward

# Basic type for differentiation.
abstract type AbstractDiff{GT, N, T} <: TagBlock{GT, N, T} end
Base.adjoint(df::AbstractDiff) = Daggered(df)
istraitkeeper(::AbstractDiff) = Val(true)
const Rotor{N, T} = Union{RotationGate{N, T}, PutBlock{N, <:Any, <:RotationGate, <:Complex{T}}}
const CphaseGate{N, T} = ControlBlock{N,<:ShiftGate{T},<:Any}
const DiffBlock{N, T} = Union{Rotor{N, T}, CphaseGate{N, T}}

# Quantum differentiation block.
"""
    QDiff{GT, N, T} <: AbstractDiff{GT, N, Complex{T}}
    QDiff(block) -> QDiff
Mark a block as quantum differentiable.
"""
mutable struct QDiff{GT, N, T} <: AbstractDiff{GT, N, Complex{T}}
    block::GT
    grad::T
    QDiff(block::DiffBlock{N, T}) where {N, T} = new{typeof(block), N, T}(block, T(0))
end

content(cb::QDiff) = cb.block
chcontent(cb::QDiff, blk::DiffBlock) = QDiff(blk)
Base.adjoint(df::QDiff) = QDiff(content(df)')
@forward QDiff.block mat, apply!

## Print the differentiation marks.
function YaoBlocks.print_annotation(io::IO, df::QDiff)
    printstyled(io, "[̂∂] "; bold=true, color=:yellow)
end

#### Functions ##### 
"""
    diff(block::AbstractBlock) -> AbstractBlock
    diff() -> Function
Mark differentiable items in a block tree(e.g.: ChainBlock) as differentiable.
"""
function diff(blk::AbstractBlock)
    blks = subblocks(blk)
    isempty(blks) ? blk : chsubblocks(blk, diff.(blks))
end
diff() = block->diff(block)
## For differentiable blocks.
diff(block::Union{RotationGate, CphaseGate}) = QDiff(block)
## Exclude control blocks.
diff(block::ControlBlock) = block



@inline function _perturb(func, gate::AbstractDiff{<:DiffBlock}, δ::Real)
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end

@inline function _perturb(func, gate::AbstractDiff{<:Rotor}, δ::Real)  # for put
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end

"""
    opdiff(psifunc, diffblock::AbstractDiff, op::MatrixBlock)
Operator differentiation.
"""
@inline function opdiff(psifunc, diffblock::AbstractDiff, op::MatrixBlock)
    r1, r2 = _perturb(()->expect(op, psifunc()) |> real, diffblock, π/2)
    diffblock.grad = (r2 - r1)/2
end