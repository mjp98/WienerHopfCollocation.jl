##--------------------------------------
## Problem = Equation + NumericalOptions
##--------------------------------------

struct ProblemWHC{N,T<:Number,A,B}
    equation::A
    numerics::B
    function ProblemWHC(equation::A, numerics::B) where {A,B}
        new{size(equation, 1),prectype(numerics),A,B}(equation, numerics)
    end
end

length(::ProblemWHC{N}) where {N} = N
prectype(::ProblemWHC{N,T}) where {N,T} = T

eltype(x::ProblemWHC) = x |> prectype |> complex
equation(x::ProblemWHC) = x.equation
numerics(x::ProblemWHC) = x.numerics

@forward ProblemWHC.equation size, evaluateF, evaluateA, evaluateB, evaluate_equation, evaluate_forcing, eachindex, eachmatrixindex
@forward ProblemWHC.numerics collocation, space, basis, nbasis, ncollocation, ifvanish, ifsave, collocationpoints, isabove, solver

# Issue with constructing PseudoBlockMatrix with tuple, so returning vector.

function blocksizes(x::ProblemWHC{N}, i::Int) where {N}
    @assert i in (1, 2)
    i == 1 && return [ncollocation(x) for _ in 1:N]
    i == 2 && return [nbasis(x) for _ in 1:N]
end

blocksizes(x::ProblemWHC) = blocksizes(x, 1), blocksizes(x, 2)

blocksize(x::ProblemWHC, i::Int) = first(blocksizes(x, i))
blocksize(x::ProblemWHC) = blocksize(x, 1), blocksize(x, 2)

function init(::Type{<:ProblemWHC}, equation; kwargs...)
    opts = OptionsWHC(; kwargs...)
    return ProblemWHC(equation, opts)
end
