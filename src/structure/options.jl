##----------------
## Options
##----------------

## Basis specification

struct BasisSpec{T}
    sp::T
    n::Int
end
nbasis(x::BasisSpec) = x.n
space(x::BasisSpec) = x.sp

@forward BasisSpec.sp isabove, prectype

## Collocation specification

struct CollocationSpec{T}
    style::T # Map
    n::Int
end
ncollocation(x::CollocationSpec) = x.n

struct OptionsWHC{B,C}
    basis::B
    collocation::C
    ifvanish::Bool # Impose vanishing condition at infinity
    ifsave::Bool # Save cauchy matrix?
    solver::Symbol
end

collocation(x::OptionsWHC) = x.collocation
basis(x::OptionsWHC) = x.basis
ifvanish(x::OptionsWHC) = x.ifvanish
ifsave(x::OptionsWHC) = x.ifsave
solver(x::OptionsWHC) = x.solver

@forward OptionsWHC.collocation ncollocation
@forward OptionsWHC.basis space, nbasis, isabove, prectype

function OptionsWHC(;
    domain=SqrtLine{-1 / 4}(0.0),
    space=Legendre(domain),
    ifvanish=true,
    ifsave=false,
    m=192,
    n=128,
    solver=:qr
)

    basis = BasisSpec(space, n)
    collocation = CollocationSpec(space, m)

    return OptionsWHC(basis, collocation, ifvanish, ifsave, solver)
end
