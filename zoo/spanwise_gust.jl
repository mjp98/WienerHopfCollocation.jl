# ------------
# Interlacing
# ------------

todoubleinterlace(n::Int) = n < 0 ? -2n : 2n + 1
fromdoubleinterlace(n::Int) = iseven(n) ? -(n ÷ 2) : (n - 1) ÷ 2 # integer division

function getinterlacedindex(x::AbstractArray{T}, i::Int) where {T}
    j = todoubleinterlace(i)
    return j <= length(x) ? x[j] : zero(T)
end

# ------------
# Spanwise varying porosity
# ------------

import WienerHopfCollocation: evaluateA, evaluateB, evaluateF

struct SpanwiseGustWH{N,T,V} <: AbstractEquationWH{N,T}
    ω::T
    M::T
    k₃::T
    compliance::V
    n::Int
end
function SpanwiseGustWH(ω, M, k₃, compliance::V, N::Int) where {V}
    args = promote(ω, M, k₃)
    T = eltype(args)
    return SpanwiseGustWH{N,T,V}(args..., compliance, N)
end

frequency(eq::SpanwiseGustWH) = eq.ω
mach_number(eq::SpanwiseGustWH) = eq.M
spanwise_wavenumber(eq::SpanwiseGustWH) = eq.k₃

function wavevector(eq::SpanwiseGustWH)
    @SVector [frequency(eq) / mach_number(eq), 1, spanwise_wavenumber(eq)]
end
wavenumber(eq::SpanwiseGustWH, i) = getindex(wavevector(eq), i)
compliance(eq::SpanwiseGustWH) = eq.compliance

shift(eqn::SpanwiseGustWH) = (size(eqn, 1) + 1) ÷ 2

function interlacedcompliance(eqn::SpanwiseGustWH, n::Int)
    return getinterlacedindex(compliance(eqn), n)
end

function evaluateA(eqn::SpanwiseGustWH, z, i, j)
    return interlacedcompliance(eqn, i - j)
end

function evaluateB(eqn::SpanwiseGustWH, z, i, j)
    k₃ = wavenumber(eqn, 3)
    k = frequency(eqn)

    w = sqrt(complex(k^2 - (k₃ + j - shift(eqn))^2))

    γ = sqrt(im * (z - w)) * sqrt(-im * (z + w))

    B = interlacedcompliance(eqn, i - j) * γ / 2
    if i == j
        B += 1
    end
    return B
end

function evaluateF(eqn::SpanwiseGustWH, z, i)
    δ = wavenumber(eqn, 1)
    k₂ = wavenumber(eqn, 2)
    if isinf(z)
        return zero(z)
    else
        return -k₂ * interlacedcompliance(eqn, i - shift(eqn)) / (z + δ)
    end
end
