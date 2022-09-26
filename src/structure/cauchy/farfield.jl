struct CauchyCache{T}
    data::Vector{T}
end

function cauchyfarfieldcache(problem,n,z=10.0)
    N = nbasis(problem)

    T = ComplexF64
    ε = 5*eps(Float64)

    d = Segment(ε,z)
    sp = space(problem)

    dm = Segment(-z,-ε)
    Cm = Matrix{T}(undef,n,N)
    unsafestieltjesmatrix!(Cm, sp, dm)

    dp = Segment(5eps(1.0),z)
    Cp = Matrix{T}(undef,n,N)
    unsafestieltjesmatrix!(Cp, sp, dp)
    d3 = Segment(-im*z,-5im*eps(1.0))
    C3 = Matrix{T}(undef,n,N)
    unsafestieltjesmatrix!(C3, sp, d3)

    d4 = Segment(5im*eps(1.0),im*z)
    C4 = Matrix{T}(undef,n,N)
    unsafestieltjesmatrix!(C4, sp, d4)

    myzero = 0.0#5eps(1.0)
    sp = space(problem)
    dm = Chebyshev(Segment(-z,-myzero))
    dp = Chebyshev(Segment(myzero,z))
    d3 = Chebyshev(Segment(-im*z,-im*myzero))
    d4 = Chebyshev(Segment(im*myzero,im*z))
    Ap = [Fun(dp,ApproxFun.transform(dp,Cp[:,i])) for i = 1:N]
    Am = [Fun(dm,ApproxFun.transform(dm,Cm[:,i])) for i = 1:N]
    A3 = [Fun(d3,ApproxFun.transform(d3,C3[:,i])) for i = 1:N]
    A4 = [Fun(d4,ApproxFun.transform(d4,C4[:,i])) for i = 1:N]
    return CauchyCache([-(Am)/(2π*im), -(Ap)/(2π*im),-(A3)/(2π*im),-(A4)/(2π*im)])
end

function unsafestieltjesmatrix!(C, sp, d)
    m = size(C,1)
    pts = points(d, m)
    for k in 1:m
        stieltjesmoment!(view(C,k,:), sp, pts[k])
    end
    C
end

function cauchycached(f::Fun, zcache)
    return cauchycached(coefficients(f), zcache)
end
function cauchycached(f::AbstractVector, zcache)
    return ComplexF64.(ApproxFun.dotu(f, zcache))
end





# function cauchyfarfieldcache(problem,m::Int,r::Real=10.0)
#     sp = space(problem)
#     d = Segment(0,r)
#     n = nbasis(problem)
#     return [unsafestieltjesmatrix(sp,cispi(i/2)*d, m, n)  for i in 1:4]
# end

# function unsafestieltjesmatrix(sp::Space, d::Segment, m::Int, n::Int)
#     T = ComplexF64
#     ε = 5*eps(real(T))

#     dom = d
#     if  abs(d.a) < 2ε
#         dom = Segment(ε*cis(angle(d.a)),d.b)
#     end
#     if abs(d.b) < 2ε
#         dom = Segment(d.a,ε*cis(angle(d.b)))
#     end

#     C = Matrix{T}(undef,m,n)
#     unsafestieltjesmatrix!(C, sp, dom)

#     A = [Fun(d,ApproxFun.transform(d,C[:,i])) for i in 1:n]

#     return -A/(2π*im)
# end
