
##--------------------------------------
## Solution wrapper
##--------------------------------------

abstract type AbstractSolutionWHC <: AbstractSolutionSplit end

problem(x:: AbstractSolutionWHC) = x.problem

@forward AbstractSolutionWHC.problem space, size, isabove

struct SolutionWHC{T,S<:Fun} <:AbstractSolutionWHC
    problem::T
    ϕ::Vector{S}
end

function SolutionWHC(problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    ϕ = [Fun(space(problem), b) for b in blocks(x)]
    return SolutionWHC(problem,ϕ)
end


function SolutionWHC!(ϕ,problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    for i in 1:blocklength(x)
        ϕ[i].coefficients .= view(x,Block(i))
    end
    return SolutionWHC(problem,ϕ)
end

function solutionfactors(solution::T) where T<: AbstractSolutionWHC
    φ₊ = SolutionFactor{true}(solution)
	φ₋ = SolutionFactor{false}(solution)
    return φ₊, φ₋
end

# Evaluation of solution

(Φ::AbstractSolutionWHC)(z...) = evaluate(Φ,z...)

cauchy(solution::T,z::Number) where T<:AbstractSolutionWHC = [cauchy(f, z) for f in solution.ϕ]

function evaluate(solution::AbstractSolutionWHC,z::T) where T<:Number
    inUHP = isabove(space(solution),z)
    cz = cauchy(solution,z)
    @unpack A,B,F = equation(solution)
    if inUHP
        return cz,  A(z) \ (F(z) - B(z) * cz)
    else
        return B(z) \ (F(z) - A(z) * cz), cz
    end
end


function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchy(solution,z)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

# If we just broadcast then get array of solutions, whereas maybe want solution of arrays?
function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:AbstractArray{<:Number}
    @unpack nwh = solution.prob
    p = [zeros(complex(float(eltype(T))),size(z)) for i = 1:nwh]
    for (i, x) in enumerate(z)
        px = solution(x,isUHP)
        for j = 1:nwh
            p[j][i] = px[j]
        end
    end
    return p
end

function evaluatecached(solution::AbstractSolutionWHC,z::T,isUHP::Bool,cache) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchycached(solution,z,cache)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

function cauchycached(solution::T,z::Complex{S},cache) where {T<:AbstractSolutionWHC,S}

    n = 1
    if (isapprox(angle(z),π)) || (isapprox(angle(z),-π))
        n = 1
    elseif angle(z) == 0
        n = 2
    elseif isapprox(angle(z),π/2)
        n = 4
    elseif isapprox(angle(z),-π/2)
        n = 3
    end

    nb = nbasis(solution.problem)
    zcache = Vector{Complex{S}}(undef,nb)
    for i = 1:nb
        zcache[i] = cache[n][i](z)
    end

    N = length(solution.ϕ)
    ret = Vector{Complex{S}}(undef,N)
    for i=1:N
        ret[i] = cauchycached(solution.ϕ[i], zcache)
    end
    return ret
end

# function cauchycached(f::T,z::S,cache) where {T<:Fun,S<:Number}
#     if z<0
#         zcache = [complex(cache[1][i](z)) for i = 1:ncoefficients(f)]
#     else
#         zcache = [complex(cache[2][i](z)) for i = 1:ncoefficients(f)]
#     end
#     return cauchycached(f::T,zcache)
#  end

function cauchycached(f::T,zcache) where {T<:Fun}
    return ComplexF64.(ApproxFun.dotu(f.coefficients,zcache))
 end

##--------------------------------------
## Solution wrapper
##--------------------------------------

abstract type AbstractSolutionWHC <: AbstractSolutionSplit end

problem(x::T) where T<: AbstractSolutionWHC = x.problem
for op in [:space,:size]
    @eval $op(x::T) where T<:AbstractSolutionWHC  = $op(x.problem)
end

isabove(solution::T,z)  where T<:AbstractSolutionWHC = isabove(space(solution),z)


struct SolutionWHC{T,S<:Fun} <:AbstractSolutionWHC
    problem::T
    ϕ::Vector{S}
end

function SolutionWHC(problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    ϕ = [Fun(space(problem), b) for b in blocks(x)]
    return SolutionWHC(problem,ϕ)
end


function SolutionWHC!(ϕ,problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    for i in 1:blocklength(x)
        ϕ[i].coefficients .= view(x,Block(i))
    end
    return SolutionWHC(problem,ϕ)
end


function solutionfactors(solution::T) where T<: AbstractSolutionWHC
    φ₊ = SolutionFactor{true}(solution)
	φ₋ = SolutionFactor{false}(solution)
    return φ₊, φ₋
end

# Evaluation of solution

(Φ::AbstractSolutionWHC)(z...) = evaluate(Φ,z...)

cauchy(solution::T,z::Number) where T<:AbstractSolutionWHC = [cauchy(f, z) for f in solution.ϕ]

function evaluate(solution::AbstractSolutionWHC,z::T) where T<:Number
    inUHP = isabove(space(solution),z)
    cz = cauchy(solution,z)
    @unpack A,B,F = equation(solution)
    if inUHP
        return cz,  A(z) \ (F(z) - B(z) * cz)
    else
        return B(z) \ (F(z) - A(z) * cz), cz
    end
end


function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchy(solution,z)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

# If we just broadcast then get array of solutions, whereas maybe want solution of arrays?
function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:AbstractArray{<:Number}
    @unpack nwh = solution.prob
    p = [zeros(complex(float(eltype(T))),size(z)) for i = 1:nwh]
    for (i, x) in enumerate(z)
        px = solution(x,isUHP)
        for j = 1:nwh
            p[j][i] = px[j]
        end
    end
    return p
end

function evaluatecached(solution::AbstractSolutionWHC,z::T,isUHP::Bool,cache) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchycached(solution,z,cache)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

function cauchycached(solution::T,z::Complex{S},cache) where {T<:AbstractSolutionWHC,S}

    n = 1
    if (isapprox(angle(z),π)) || (isapprox(angle(z),-π))
        n = 1
    elseif angle(z) == 0
        n = 2
    elseif isapprox(angle(z),π/2)
        n = 4
    elseif isapprox(angle(z),-π/2)
        n = 3
    end

    nb = nbasis(solution.problem)
    zcache = Vector{Complex{S}}(undef,nb)
    for i = 1:nb
        zcache[i] = cache[n][i](z)
    end

    N = length(solution.ϕ)
    ret = Vector{Complex{S}}(undef,N)
    for i=1:N
        ret[i] = cauchycached(solution.ϕ[i], zcache)
    end
    return ret
end

# function cauchycached(f::T,z::S,cache) where {T<:Fun,S<:Number}
#     if z<0
#         zcache = [complex(cache[1][i](z)) for i = 1:ncoefficients(f)]
#     else
#         zcache = [complex(cache[2][i](z)) for i = 1:ncoefficients(f)]
#     end
#     return cauchycached(f::T,zcache)
#  end

function cauchycached(f::T,zcache) where {T<:Fun}
    return ComplexF64.(ApproxFun.dotu(f.coefficients,zcache))
 end
##--------------------------------------
## Solution wrapper
##--------------------------------------

abstract type AbstractSolutionWHC <: AbstractSolutionSplit end

problem(x::T) where T<: AbstractSolutionWHC = x.problem
for op in [:space,:size]
    @eval $op(x::T) where T<:AbstractSolutionWHC  = $op(x.problem)
end

isabove(solution::T,z)  where T<:AbstractSolutionWHC = isabove(space(solution),z)


struct SolutionWHC{T,S<:Fun} <:AbstractSolutionWHC
    problem::T
    ϕ::Vector{S}
end

function SolutionWHC(problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    ϕ = [Fun(space(problem), b) for b in blocks(x)]
    return SolutionWHC(problem,ϕ)
end


function SolutionWHC!(ϕ,problem::ProblemWHC,x::PseudoBlockVector{T}) where T<:Number
    for i in 1:blocklength(x)
        ϕ[i].coefficients .= view(x,Block(i))
    end
    return SolutionWHC(problem,ϕ)
end


function solutionfactors(solution::T) where T<: AbstractSolutionWHC
    φ₊ = SolutionFactor{true}(solution)
	φ₋ = SolutionFactor{false}(solution)
    return φ₊, φ₋
end

# Evaluation of solution

(Φ::AbstractSolutionWHC)(z...) = evaluate(Φ,z...)

cauchy(solution::T,z::Number) where T<:AbstractSolutionWHC = [cauchy(f, z) for f in solution.ϕ]

function evaluate(solution::AbstractSolutionWHC,z::T) where T<:Number
    inUHP = isabove(space(solution),z)
    cz = cauchy(solution,z)
    @unpack A,B,F = equation(solution)
    if inUHP
        return cz,  A(z) \ (F(z) - B(z) * cz)
    else
        return B(z) \ (F(z) - A(z) * cz), cz
    end
end


function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchy(solution,z)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

# If we just broadcast then get array of solutions, whereas maybe want solution of arrays?
function evaluate(solution::AbstractSolutionWHC,z::T,isUHP::Bool) where T<:AbstractArray{<:Number}
    @unpack nwh = solution.prob
    p = [zeros(complex(float(eltype(T))),size(z)) for i = 1:nwh]
    for (i, x) in enumerate(z)
        px = solution(x,isUHP)
        for j = 1:nwh
            p[j][i] = px[j]
        end
    end
    return p
end

function evaluatecached(solution::AbstractSolutionWHC,z::T,isUHP::Bool,cache) where T<:Number
    @unpack A, B, F = equation(solution.problem)
    inUHP = isabove(solution,z)
    cz = cauchycached(solution,z,cache)
    if isUHP
        return inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
    else
        return inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
    end
end

function cauchycached(solution::T,z::Complex{S},cache) where {T<:AbstractSolutionWHC,S}

    n = 1
    if (isapprox(angle(z),π)) || (isapprox(angle(z),-π))
        n = 1
    elseif angle(z) == 0
        n = 2
    elseif isapprox(angle(z),π/2)
        n = 4
    elseif isapprox(angle(z),-π/2)
        n = 3
    end

    nb = nbasis(solution.problem)
    zcache = Vector{Complex{S}}(undef,nb)
    for i = 1:nb
        zcache[i] = cache[n][i](z)
    end

    N = length(solution.ϕ)
    ret = Vector{Complex{S}}(undef,N)
    for i=1:N
        ret[i] = cauchycached(solution.ϕ[i], zcache)
    end
    return ret
end

# function cauchycached(f::T,z::S,cache) where {T<:Fun,S<:Number}
#     if z<0
#         zcache = [complex(cache[1][i](z)) for i = 1:ncoefficients(f)]
#     else
#         zcache = [complex(cache[2][i](z)) for i = 1:ncoefficients(f)]
#     end
#     return cauchycached(f::T,zcache)
#  end

function cauchycached(f::T,zcache) where {T<:Fun}
    return ComplexF64.(ApproxFun.dotu(f.coefficients,zcache))
 end
