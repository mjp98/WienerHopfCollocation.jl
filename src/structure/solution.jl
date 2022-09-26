
##--------------------------------------
## Solution wrapper
##--------------------------------------

abstract type AbstractSplit <: Function end
abstract type AbstractSolutionSplit <: AbstractSplit end

abstract type AbstractSolutionWH{T} <: AbstractSolutionSplit end

problem(x::AbstractSolutionWH) = x.problem
@forward AbstractSolutionWH.problem space, size, isabove
@forward AbstractSolutionWH.problem evaluateA, evaluateB, evaluateF
@forward AbstractSolutionWH.problem evaluate_equation, evaluate_forcing

(Φ::AbstractSolutionWH)(z...) = evaluate(Φ, z...)

struct SolutionWH{T<:Number,R,S} <: AbstractSolutionWH{T}
    problem::R
    density::Vector{S}
    function SolutionWH(problem::A, density::Vector{B}) where {A,B}
        return new{ComplexF64,A,B}(problem, density)
    end
end

function SolutionWH!(phi, problem::ProblemWHC, x::PseudoBlockVector{<:Number})
    for i in 1:blocklength(x)
        phi[i].coefficients .= view(x, Block(i))
    end
    return SolutionWH(problem, phi)
end

SolutionWH!(sol, problem::ProblemWHC) = SolutionWH!(sol.phi, problem, sol.x)

density(x::SolutionWH) = x.density
coefficients(x::SolutionWH, i) = coefficients(getindex(density(x), i))
cauchy(x::AbstractSolutionWH, z) = [cauchy(f, z) for f in density(x)]
length(x::SolutionWH) = x |> density |> length

function evaluate(x::AbstractSolutionWH, z::Number, C)
    A, B = evaluate_equation(x, z)
    F = evaluate_forcing(x, z)

    if isabove(x, z)
        return C, B \ (F - A * C)
    else
        return A \ (F - B * C), C
    end
end

function evaluate(x::AbstractSolutionWH, z, c::CauchyCache)
    return evaluate(x, z, cauchy(x, z, c))
end
evaluate(x::AbstractSolutionWH, z) = evaluate(x, z, cauchy(x, z))


function evaluate(x::AbstractSolutionWH, z::Number, u::Bool, C)
    A, B = evaluate_equation(x, z)
    F = evaluate_forcing(x, z)

    inupperhalfplane = isabove(x, z)

    (u  &&  inupperhalfplane) && return C
    (u  && !inupperhalfplane) && return B\(F-A*C)
    (!u && !inupperhalfplane) && return C
    (!u &&  inupperhalfplane) && return A\(F-B*C)
end

function evaluate(x::AbstractSolutionWH, z::Number, u::Bool, c::CauchyCache)
    return evaluate(x, z, u, cauchy(x, z, c))
end
function evaluate(x::AbstractSolutionWH, z::Number, u::Bool)
    return evaluate(x, z, u, cauchy(x, z))
end

function cauchy(solution::AbstractSolutionWH, z, cache::CauchyCache)

    T = ComplexF64
    if isapprox(abs(angle(z)), π)
        n = 1
    elseif isapprox(angle(z), 0)
        n = 2
    elseif isapprox(angle(z), -π / 2)
        n = 3
    elseif isapprox(angle(z), π / 2)
        n = 4
    else
        @warn "cauchy(solution, z) not in cache"
        return cauchy(solution, z)
    end

    nb = nbasis(solution.problem)
    zcache = Vector{T}(undef, nb)
    for i in 1:nb
        zcache[i] = cache.data[n][i](z)
    end

    N = length(solution)
    ret = Vector{T}(undef, N)
    for i in 1:N
        ret[i] = cauchycached(coefficients(solution, i), zcache)
    end
    return ret
end
