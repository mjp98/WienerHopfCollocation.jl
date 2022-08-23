module WienerHopfCollocation

using Reexport

@reexport using ApproxFun
using BlockArrays
@reexport using CommonSolve
using FastLapackInterface
using JLD2
using Lazy
using LinearAlgebra
using RiemannHilbert
@reexport using ScaledLines
@reexport using SingularIntegralEquations
using StaticArrays
using UnPack
@reexport using WienerHopf


import Base:size
import RiemannHilbert: collocationpoints

export isabove, preallocate, precompute
export WienerHopfEq
export collocationrhs, collocationmatrix, collocationmatrix!, collocationrhs!

include("problem.jl")
include("split.jl")
include("solution.jl")
include("cauchy.jl")

collocationpoints(x::ProblemWHC) = collocationpoints(space(x),collocation(x))
collocationpoints(x::Space,s::CollocationSpec) = points(x,ncollocation(s))


# Rewrite

# TODO: organise this code
# TODO: package this code
# TODO: tests and Coverage.jl this code
# TODO: document Documenter.jl this code
# TODO: Profile.jl and ProfileViews.jl this code
# TODO: Plot Line() Tikzpicture, PlotRay()
# TODO: Want structure that include problem parameters and solution...

struct HelmholtzKernel{T<:Complex}
    k::T
end

function γm(x,k)
    return sqrt(im*(x-k))
end
function γp(x,k)
    return sqrt(-im*(x+k))
end
function helmholtzγ(x,k)
    return  γm(x,k)*γp(x,k)
end


cfeltype(d::Domain) = complex(float(eltype(d)))
cfeltype(S::Space) = cfeltype(domain(S))

## isabove

isabove(sp::Space,z::T) where T<:Number = isabove(domain(sp), z)
isaboveline(c,θ,z) = imag(cis(-θ)*(z - c)) > 0
isabove(dom::Union{Line,SqrtLine},z) = isaboveline(dom.center,angle(dom),z)
isabove(dom::ScaledSegmentLine,z) = isabove(descale(dom),z)


##--------------------------------------
## Cauchy Matrix
##--------------------------------------

function fpcauchymatrix2(x...)
    C = fpstieltjesmatrix2(x...)
    C./=π
    C ./= (-2*im)
    C
end

function fpstieltjesmatrix2(spec::BasisSpec,pts::Vector{T}) where T
    @unpack n,sp = spec
    return fpstieltjesmatrix2(sp, n, pts)
end
function fpstieltjesmatrix2(sp::Space, n::Int, pts::Vector{T}) where T
    m = length(pts)
    return fpstieltjesmatrix2!(Array{T}(undef, m, n), sp, domain(sp),pts)
end

# Use symmetry
function fpstieltjesmatrix2!(C, sp, d, pts2::Vector{T}) where T
    pts = deepcopy(pts2)
    m, n = size(C)
    if d == domain(sp)
        #fprightstieltjesmoment!(view(C,1,:), sp)
        for k=2:Int(ceil(m/2))
            stieltjesmoment!(view(C,k,:), sp, Directed{false}(pts[k]))
        end
       # @info "assuming symmetry:\n C[m+1-k,2j-1] = -conj(C[k,2j-1]) and C[m+1-k,2j] = conj(C[k,2j])"

        # TODO: change this to a single loop.
        for k = 1:Int(floor(m/2))
            for j = 2:2:n
                C[m+1-k,j] = conj(C[k,j])
            end
            for j = 1:2:n
                C[m+1-k,j] = -conj(C[k,j])
            end
        end
    elseif leftendpoint(d) ∈ domain(sp) && rightendpoint(d) ∈ domain(sp)
        fprightstieltjesmoment!(view(C,1,:), sp, d)
        for k=2:m-1
            stieltjesmoment!(view(C,k,:), sp, pts[k])
        end
        fpleftstieltjesmoment!(view(C,m,:), sp, d)
    elseif leftendpoint(d) ∈ domain(sp)
        for k=1:m-1
            stieltjesmoment!(view(C,k,:), sp, pts[k])
        end
        fpleftstieltjesmoment!(view(C,m,:), sp, d)
    elseif rightendpoint(d) ∈ domain(sp)
        fprightstieltjesmoment!(view(C,1,:), sp, d)
        for k=2:m
            stieltjesmoment!(view(C,k,:), sp, pts[k])
        end
    else
        for k=1:m
            stieltjesmoment!(view(C,k,:), sp, pts[k])
        end
    end
    C
end

function myevaluationmatrix(basis::BasisSpec,pts::Vector{T}) where T
    @unpack n,sp = basis
    m = length(pts)
    E = Matrix{T}(undef,m,n)
    return evaluationmatrix!(E, sp,deepcopy(pts))
end

function filename_cauchymatrix(sp,m,n)
    T = string(real(cfeltype(sp)))
    if real(cfeltype(sp)) <: BigFloat
        @warn "precision of bigfloats not saved - if loading old data and have changed precision, delete file to recompute"
    end
    sc = string(Int(log10(scale(sp))))
    !isdirpath("cauchydata") && mkdir("cauchydata")
    filename = "cauchydata/cauchymatrix_m"*string(Int(m))*"_n"*string(Int(n))*"_oftype_"*T*"scale_"*sc*".jld2"
    return filename
end

function construct_cauchymatrix(problem::ProblemWHC,pts)
   # C₋ = fpcauchymatrix2(basis(problem),pts)
    C₋ = fpcauchymatrix(space(problem),domain(space(problem)),ncollocation(problem),nbasis(problem))
    E  = myevaluationmatrix(basis(problem),pts)
    if ifvanish(problem)
        E  = ensure_vanish!(E)
        C₋ = ensure_vanish!(C₋)
    end
    C₊ = C₋ + E
    return C₋,C₊
end

struct CauchyData{T}
    C₋::Matrix{T}
    C₊::Matrix{T}
    pts::Vector{T}
end

function makecauchy(problem,pts;ifverbose=false)
    if ifsave(problem)
        filename = filename_cauchymatrix(space(problem),length(pts),nbasis(problem))
        if isfile(filename)
            ifverbose && @info "loading "*filename
            C₋,C₊= load(filename,"Cm", "Cp")
        else
            ifverbose && @info "no data found at "*filename
            C₋,C₊ = construct_cauchymatrix(problem,pts)

            ifverbose && @info "saving "*filename
            save(filename, "Cm",C₋,"Cp",C₊)
        end
    else
        C₋,C₊ = construct_cauchymatrix(problem,pts)
    end
    return CauchyData(C₋,C₊,pts)
end

	# precompute


	function precompute(problem)
        pts = collocationpoints(problem)
        return makecauchy(problem,pts)
    end

##--------------------------------------
## Preallocation
##--------------------------------------

# TODO: neaten up how type is inferred... (currently from pts...)

function preallocate_matrix(problem::ProblemWHC,::Vector{T},mwh,nwh) where T
    m = @SVector [ncollocation(problem)]
    n = @SVector [nbasis(problem)]
    return PseudoBlockMatrix{T}(undef,repeat(m, mwh), repeat(n, nwh))
end

function preallocate_rhs(problem::ProblemWHC,::Vector{T},mwh) where T
    m = @SVector [ncollocation(problem)]
    return PseudoBlockVector{T}(undef, repeat(m, mwh))
end

function preallocate_lhs(problem::ProblemWHC,::Vector{T},nwh) where T
    n = @SVector [nbasis(problem)]
    return PseudoBlockVector{T}(undef, repeat(n, nwh))
end

function preallocate_sol(problem::ProblemWHC,::Vector{T}) where T
	N,_ = size(problem)
	n  = nbasis(problem)
	sp = space(problem)
	F = [Fun(sp,zeros(T,n)) for i = 1:N]
	return F
end

function preallocate(problem::ProblemWHC)
    pts = collocationpoints(problem)
    (mwh,nwh) = size(problem)
    return preallocate(problem,pts,mwh,nwh)
end

function preallocate(problem::ProblemWHC,pts::Vector{T},mwh,nwh) where T
    x = preallocate_lhs(problem,pts,nwh)
    A = preallocate_matrix(problem,pts,mwh,nwh)
    b = preallocate_rhs(problem,pts,mwh)
    φ = preallocate_sol(problem,pts)
    return WHCdata(x,A,b,φ)
end



##--------------------------------------
## Collocation system
##--------------------------------------

function collocationmatrix(problem::ProblemWHC,pts)
    cdata = makecauchy(problem,pts)
    return collocationmatrix(problem,pts,cdata)
end

function collocationmatrix(problem::ProblemWHC,pts,cdata)
    ret = preallocate_matrix(problem,pts)
    return collocationmatrix!(ret,problem,pts,cdata)
end

function collocationmatrix!(ret,problem::ProblemWHC,pts::Vector{T},C::CauchyData) where T
    (_,nwh) = size(problem)
    return pseudocollocationmatrixvanish!(ret,pts,nwh,
                    problem.equation.A,C.C₋,
                    problem.equation.B,C.C₊)
end

function collocationmatrix!(ret,problem::ProblemWHC,C::CauchyData) where T
    (_,nwh) = size(problem)
    return pseudocollocationmatrixvanish!(ret,C.pts,nwh,
                    problem.equation.A,C.C₋,
                    problem.equation.B,C.C₊)
end

function pseudocollocationmatrix!(ret,pts::Vector{T},nwh,A,C₋,B,C₊) where T
    m, n = size(C₊)
    fill!(ret,zero(T))
    for J = 1:nwh
        for I = 1:nwh
            for i=1:m
                x = pts[i]
                a = A(x,I,J)
                b = B(x,I,J)
                for j = 1:n
                    ret[BlockIndex((I,J), (i,j))] = a*C₋[i,j] + b*C₊[i,j]
                end
            end
        end
    end
    return ret
end

function pseudocollocationmatrixvanish!(ret,pts::Vector{T},nwh,A,C₋,B,C₊) where T
    # TODO: better exploit column ordering
    m, n = size(C₊)
    fill!(ret,zero(eltype(ret)))
    for J = 1:nwh
        for I = 1:nwh
            # first row
            if I==J
                for j = 1:n
                    ret[BlockIndex((I,J), (1,j))] = one(T)
                end
            end
            for i=2:(m-1)
                x = pts[i]
                a = A(x,I,J)
                b = B(x,I,J)
                for j = 1:n
                    ret[BlockIndex((I,J), (i,j))] = a*C₋[i,j] + b*C₊[i,j]
                end
            end
            # last row
            if I==J
                for j = 1:n
                    ret[BlockIndex((I,J), (m,j))] = one(T)*sign(mod(j,2)-1/2)
                end
            end
        end
    end
    return ret
end

##--------------------------------------
## Collocation RHS
##--------------------------------------

function collocationrhs(problem::ProblemWHC,pts)
    b = preallocate_rhs(problem,pts)
    return collocationrhs!(b,problem::ProblemWHC,pts)
end

function collocationrhs!(b,problem::ProblemWHC,cdata::CauchyData)
    (_,nwh) = size(problem)
    return forcingvanish!(b,cdata.pts,nwh,problem.equation.F)
end

function collocationrhs!(b,problem::ProblemWHC,pts)
    (_,nwh) = size(problem)
    return forcingvanish!(b,pts,nwh,problem.equation.F)
end

function forcing!(ret,pts::Vector{T},nwh,F) where {T<:Complex}
    n = length(pts)
    for i = 1:nwh
        for k = 1:n
            ret[BlockIndex((i),(k))] = F(pts[k],i)
        end
    end
    return ret
end

function forcingvanish!(ret,pts::Vector{T},nwh,F) where {T<:Complex}
    n = length(pts)
    for i = 1:nwh
        ret[BlockIndex((i),(1))] = zero(T)
        for k = 2:(n-1)
            ret[BlockIndex((i),(k))] = F(pts[k],i)
        end
        ret[BlockIndex((i),(n))] = zero(T)
    end
    return ret
end

##--------------------------------------
## Problem wrapper
##--------------------------------------

function SolveSpec(;domain=SqrtLine{-1/4}(zero(ComplexF64)),space=:Legendre,ifvanish=true,ifsave=true,m=80,n=80,solver=:lu)
    expansionspace = BasisSpec(space,n)
    collocationoptions = CollocationSpec(:space,m)
    return SolveSpec(expansionspace,collocationoptions,ifvanish,ifsave,solver)
end

function makeproblem(equation::WienerHopfEq,d=SqrtLine{-1/4}(zero(ComplexF64));kwargs...)
    numopts = numericaloptions(d;kwargs...)
    return ProblemWHC(equation,numericaloptions)
end


function problem(eqn,numerical_options)
    opts = SolveSpec(;numerical_options...)
    return ProblemWHC(eqn,opts)
end


##--------------------------------------
## Solve wrapper
##--------------------------------------

# TODO: put pts into cdata
# TODO: somehow make Vector{Fun} use Blockarray for memory

function solve(spec::ProblemWHC)
    preallocated = preallocate(spec);
    precomputed = precompute(spec);
    return lusolve!(preallocated,precomputed,spec)
end

# solve(problem::ProblemWHC) = solve(problem,collocationpoints(problem))
# solve(problem::ProblemWHC,pts) = solve(problem,pts,makecauchy(problem,pts))

# function solve(problem::ProblemWHC,cdata::CauchyData)

#     sdata = preallocate(problem)

# 	return solve!(sdata,problem,cdata)
# end

function solve!(sdata::WHCdata,problem::ProblemWHC,cdata)
    pts = cdata.pts
    A = collocationmatrix!(sdata.A,problem,pts,cdata)
    b = collocationrhs!(sdata.b,problem,pts)
	return solve!(sdata,problem)
end

function solve!(sdata::WHCdata,cdata,problem::ProblemWHC)
    if solver(problem) == :lu
        return lusolve!(sdata,cdata,problem)
    elseif solver(problem) == :recursivelu
        return recursivelusolve!(sdata,cdata,problem)
    elseif sovler(problem) == :sparse
        return sparsesolve!(sdata,cdata,problem)
    else # solver(problem) == :qr
        return qrsolve!(sdata,cdata,problem)
    end
end

function qrsolve!(sdata::WHCdata,cdata,problem::ProblemWHC)
    pts = cdata.pts
    collocationmatrix!(sdata.A,problem,pts,cdata)
    collocationrhs!(sdata.b,problem,pts)
    qrsolve!(sdata)
	return SolutionWHC!(sdata.φ,problem,sdata.x)
end

function lusolve!(sdata::WHCdata,cdata,problem::ProblemWHC)
    pts = cdata.pts
    collocationmatrix!(sdata.A,problem,pts,cdata)
    collocationrhs!(sdata.b,problem,pts)
    lusolve!(sdata)
	return SolutionWHC!(sdata.φ,problem,sdata.x)
end

# function recursivelusolve!(sdata::WHCdata,cdata,problem::ProblemWHC)
#     pts = cdata.pts
#     collocationmatrix!(sdata.A,problem,pts,cdata)
#     collocationrhs!(sdata.b,problem,pts)
#     recursivelusolve!(sdata)
# 	return SolutionWHC!(sdata.φ,problem,sdata.x)
# end

function sparsesolve!(sdata::WHCdata,cdata,problem::ProblemWHC)
    pts = cdata.pts
    collocationmatrix!(sdata.A,problem,pts,cdata)
    collocationrhs!(sdata.b,problem,pts)
    sparsesolve!(sdata)
	return SolutionWHC!(sdata.φ,problem,sdata.x)
end

function qrsolve!(sdata::WHCdata)
    temp = ldiv!(qr!(sdata.A.blocks),sdata.b.blocks)
    sdata.x.blocks .=  @views temp[1:length(sdata.x)]
end

function lusolve!(sdata::WHCdata)
    temp = ldiv!(lu!(sdata.A.blocks),sdata.b.blocks)
    sdata.x.blocks .=  @views temp[1:length(sdata.x)]
end

function sparsesolve!(sdata::WHCdata)
    sdata.x.blocks .= sparse(sdata.A.blocks)\sdata.b.blocks
end

# function recursivelusolve!(sdata::WHCdata)
#     sdata.x.blocks .= ldiv!(RecursiveFactorization.lu!(sdata.A.blocks),sdata.b.blocks)
# end

# JacobiOrPolynomialSqrtLine = Union{PolynomialSpace{<:SqrtLine},JacobiWeight{<:PolynomialSpace{<:SqrtLine}}}


# function stieltjesmoment!(ret,S::JacobiOrPolynomialSqrtLine,z,filter=identity)
#     if domain(S) == SqrtLine()
#         s = setcanonicaldomain(S)
#         u = complex(undirected(z))
#         tmp = similar(ret)
#         stieltjesmoment!(ret, s, ifourtoone(1,z),filter)
#         stieltjesmoment!(tmp, s, ifourtoone(2,u),filter); ret .+= tmp
#         stieltjesmoment!(tmp, s, ifourtoone(3,u),filter); ret .+= tmp
#         stieltjesmoment!(tmp, s, ifourtoone(4,u),filter); ret .+= tmp
#         fpleftstieltjesmoment!(tmp, s); ret .-= 2 .* tmp
#         fprightstieltjesmoment!(tmp, s); ret .-= 2 .* tmp
#         ret
#     else
#         stieltjesmoment!(ret,setdomain(S,SqrtLine()),mappoint(domain(S),SqrtLine(),z),filter)
#     end
# end


# fprightstieltjesmoment!(V, sp::JacobiOrPolynomialSqrtLine) = fill!(V,0)
# fpleftstieltjesmoment!(V, sp::JacobiOrPolynomialSqrtLine) = fill!(V,0)


# import SingularIntegralEquations:directed_₂F₁general,_₂F₁maclaurin,_₂F₁Inf,_₂F₁one

# function directed_₂F₁general(a::Number,b::Number,c::Number,z)
#     T = promote_type(typeof(a),typeof(b),typeof(c),typeof(undirected(z)))

#     real(b) < real(a) && (return _₂F₁general(b,a,c,z))
#     real(c) < real(a)+real(b) && (return exp((c-a-b)*log1p(-z))*_₂F₁general(c-a,c-b,c,z))

#     #ρ = max(p, q)+1
#     ρ = max(2, 1)+1

#     if abs(z) ≤ ρ || -a ∈ ℕ₀ || -b ∈ ℕ₀
#         _₂F₁maclaurin(a,b,c,undirected(z))
#     elseif abs(z/(z-1)) ≤ ρ
#         exp(-a*log1p(-z))_₂F₁maclaurin(a,c-b,c,undirected(z/(z-1)))
#     elseif abs(inv(z)) ≤ ρ
#         _₂F₁Inf(a,b,c,z)
#     elseif abs(1-inv(z)) ≤ ρ
#         exp(-a*log1p(-z))*_₂F₁Inf(a,c-b,c,reverseorientation(z/(z-1)))
#     elseif abs(1-z) ≤ ρ
#         _₂F₁one(a,b,c,z)
#     elseif abs(inv(1-z)) ≤ ρ
#         exp(-a*log1p(-z))*_₂F₁one(a,c-b,c,reverseorientation(z/(z-1)))
#     else
#         _₂F₁taylor(a,b,c,undirected(z))
#     end
# end

end
