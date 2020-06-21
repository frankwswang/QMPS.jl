export markDiff, getQdiff!, getNdiff!
import StatsBase: mean


#### Structures #####
# Define supertypes of differentiable blocks.
const CphaseGate{N, T} = ControlBlock{N, <:ShiftGate{T}, <:Any}
const DiffBlock{N, T} = Union{RotationGate{N, T, <:Any}, CphaseGate{N, T}}


# Define the struct of a differentiable gate `QDiff{GT, N}`.
"""
    QDiff{GT, N} <: TagBlock{GT, N}
Differentiable gate. If you want to mark or return a differentiable block, use `markDiff` or `markDiff!` instead.
\n
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
    getQdiff!(psifunc::Function, diffGate::Union{QMPS.QDiff, Array{QMPS.QDiff,1}}, op::AbstractBlock) -> grad::Union{Float64, Array{Float64, 1}}
Quantum Operator differentiation. When only apply `getQdiff!` to one differentiable gate, the output satisfies `grad == diffGate.grad`.  
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `diffGate = collect_blocks(QMPS.QDiff, c)`
\n `op`: Witness Operator to measure reg.
\n 
\nNOTE: 
\n1) `getQdiff!` modifies only the element `grad` inside the field of `QDiff` and all other elements remain the same.  
\n2) `getQdiff!` does indirectly modify the input qubit `reg` to `reg |> c` when inputting the `psifunc` as an arguement. 
"""
@inline function getQdiff!(psifunc::Function, diffGate::QDiff, op::AbstractBlock)
    r1, r2 = _perturb( ()->mean( Yao.expect(op, psifunc()) ) |> real, diffGate, π/2 )
    diffGate.grad = (r2 - r1)/2
end
@inline function getQdiff!(psifunc::Function, diffGates::Array{QDiff,1}, op::AbstractBlock)
    grads = getQdiff!.(psifunc, diffGates, Ref(op))
end


"""
    getNdiff!(psifunc::Function, parGate::Union{AbstractBlock, Array{AbstractBlock,1}}, op::AbstractBlock; δ::Real=0.01) -> grad::Union{Float64, Array{Float64, 1}}
Numerical Operator differentiation. When only apply `getNdiff!!` to one differentiable gate, the output satisfies `grad == parGate.grad`.
\n `psifunc = ()-> reg::ArrayReg |> c::ChainBlock`
\n `parGate`: Paramterized block(gate) in c.
\n `op`: Witness Operator to measure reg.
\n
\nNOTE: 
\n1) When applying `getNdiff!` to QMPS-DC or similar differentiable MPS circuits, the accuracy of the function can get much worse if you set `δ` to be too small (<1e-3).
\n2) `getNdiff!` indirectly modifies the input qubit `reg` to `reg |> c` when inputting the `psifunc` as an arguement.
"""
@inline function getNdiff!(psifunc::Function, parGate::AbstractBlock, op::AbstractBlock; δ::Real=0.01)
    r1, r2 = _perturb( ()->mean( Yao.expect(op, psifunc()) ) |> real, parGate, δ )
    grad = (r2 - r1) / (2δ)
end
@inline function getNdiff!(psifunc::Function, parGates::Array{<:AbstractBlock,1}, op::AbstractBlock; δ::Real=0.01)
    grads = getNdiff!.(psifunc, parGates, Ref(op), δ=δ)
end

@inline function _perturb(func, gate::AbstractBlock, δ::Real)
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end

