##----------------
## Data structures
##----------------

struct WHCdata{T<:Number,S}
    x::PseudoBlockVector{T}
    A::PseudoBlockMatrix{T}
    b::PseudoBlockVector{T}
    φ::Vector{S}
end

##----------------
## Equation
##----------------

struct WienerHopfEq{N,M}
    A::Function
    B::Function
    F::Function
    function WienerHopfEq(A,B,F)
        (n,m) = size(A(0))
        return new{n,m}(A,B,F)
    end
end
size(::WienerHopfEq{n,m}) where {n,m} = (n,m)

##----------------
## Options
##----------------

## Basis specification

struct BasisSpec{A}
    sp::A
    n::Int
end
nbasis(x::BasisSpec) = x.n
space(x::BasisSpec) = x.sp

## Collocation specification

struct CollocationSpec{A}
    style::A # Map
    n::Int
end
ncollocation(x::CollocationSpec) = x.n


struct SolveSpec{B,C}
    basis::B
    collocation::C
    ifvanish::Bool # Impose vanishing condition at infinity
    ifsave::Bool # Save cauchy matrix?
    solver::Symbol
end

collocation(x::SolveSpec) = x.collocation
basis(x::SolveSpec) = x.basis
ifvanish(x::SolveSpec) = x.ifvanish
ifsave(x::SolveSpec) = x.ifsave
solver(x::SolveSpec) = x.solver

@forward SolveSpec.collocation ncollocation
@forward SolveSpec.basis space, nbasis

##--------------------------------------
## Problem = Equation + NumericalOptions
##--------------------------------------

struct ProblemWHC{A,B}
    equation::A
    numerics::B
end

equation(x::ProblemWHC) = x.equation
numerics(x::ProblemWHC) = x.numerics

@forward ProblemWHC.equation Base.size
@forward ProblemWHC.numerics collocation, space, basis, nbasis, ncollocation, ifvanish, ifsave
