##--------------------------------------
## Preallocation
##--------------------------------------

mutable struct DataWHC{T<:Number,S<:Fun,R}
    x::PseudoBlockVector{T}
    A::PseudoBlockMatrix{T}
    b::PseudoBlockVector{T}
    phi::Vector{S}
    ws::R
end

function preallocate(problem::ProblemWHC)
    A = preallocate_matrix(problem)
    x = preallocate_lhs(problem)
    b = preallocate_rhs(problem)
    phi = preallocate_sol(problem)
    ws = preallocate_ws(problem, A.blocks)
    return DataWHC(x, A, b, phi, ws)
end

function preallocate_matrix(x::ProblemWHC)
    return PseudoBlockMatrix{eltype(x)}(undef, blocksizes(x, 1), blocksizes(x, 2))
end
function preallocate_rhs(x::ProblemWHC)
    return PseudoBlockVector{eltype(x)}(undef, blocksizes(x, 1))
end
function preallocate_lhs(x::ProblemWHC)
    return PseudoBlockVector{eltype(x)}(undef, blocksizes(x, 2))
end
function preallocate_sol(x::ProblemWHC)
    return [Fun(space(x), zeros(eltype(x), nbasis(x))) for _ in 1:size(x, 1)]
end
function preallocate_ws(x::ProblemWHC, A)
    isbig(prectype(x)) && return nothing
    solver(x) == :qr && return Workspace(LAPACK.geqrf!, A)
    solver(x) == :lu && return Workspace(LAPACK.getrf!, A)
    @error "not implemented"
end
