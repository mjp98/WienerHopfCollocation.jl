function cauchyfarfieldcache(problem,n,z=10.0)
    N = nbasis(problem)

    sp = space(problem)
    dm = Segment(-z,-5eps(1.0))
    Cm = Matrix{ComplexF64}(undef,n,N)
    unsafestieltjesmatrix!(Cm, sp, dm)

    dp = Segment(5eps(1.0),z)
    Cp = Matrix{ComplexF64}(undef,n,N)
    unsafestieltjesmatrix!(Cp, sp, dp)


    sp = space(problem)
    d3 = Segment(-im*z,-5im*eps(1.0))
    C3 = Matrix{ComplexF64}(undef,n,N)
    unsafestieltjesmatrix!(C3, sp, d3)

    d4 = Segment(5im*eps(1.0),im*z)
    C4 = Matrix{ComplexF64}(undef,n,N)
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
    return [-(Am)/(2π*im), -(Ap)/(2π*im),-(A3)/(2π*im),-(A4)/(2π*im)]
end
