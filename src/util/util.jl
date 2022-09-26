isbig(x) = false
isbig(x::Number) = real(x) isa Union{BigFloat, BigInt}
isbig(x::Type{<:Number}) = real(x) <: Union{BigFloat, BigInt}

cfeltype(d) = complex(float(eltype(d)))
cfeltype(S::Space) = cfeltype(domain(S))


## isabove

isabove(sp::Space, z) = isabove(domain(sp), z)
isaboveline(c, θ, z) = imag(cis(-θ) * (z - c)) > 0
isabove(dom::Union{Line,SqrtLine}, z) = isaboveline(dom.center, angle(dom), z)
isabove(dom::ScaledSegmentLine, z) = isabove(ScaledLines.descale(dom), z)

# Gamma kernel

struct HelmholtzKernel{T<:Complex}
    k::T
end

γm(x, k) = sqrt( im * (x - k))
γp(x, k) = sqrt(-im * (x + k))
helmholtzγ(x, k) = γm(x, k) * γp(x, k)
