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

struct SpanwiseWH{N,T,M,V} <: AbstractEquationWH{N,T}
    incident::PlaneWave{M,T}
    compliance::V
    n::Int
end
function SpanwiseWH(wave::PlaneWave{M,T}, compliance::V, N::Int) where {M,T,V}
    return SpanwiseWH{N,T,M,V}(wave, compliance, N)
end

function SpanwiseWH{N}(wave::PlaneWave{M,T}, compliance::V) where {N,M,T,V}
    return SpanwiseWH{N,T,M,V}(wave, compliance, N)
end
function SpanwiseWH(wave::PlaneWave{M,T}, compliance::V) where {M,T,V}
    return SpanwiseWH{length(compliance)}(wave, compliance)
end

incident(eq::SpanwiseWH) = eq.incident
frequency(eq::SpanwiseWH) = eq |> incident |> frequency
wavevector(eq::SpanwiseWH) = eq |> incident |> wavevector
wavenumber(eq::SpanwiseWH, i) = getindex(wavevector(eq), i)
compliance(eq::SpanwiseWH) = eq.compliance

shift(eqn::SpanwiseWH) = (size(eqn, 1) + 1) ÷ 2

function interlacedcompliance(eqn::SpanwiseWH, n::Int)
    return getinterlacedindex(compliance(eqn), n)
end

function evaluateA(eqn::SpanwiseWH, z, i, j)
    return interlacedcompliance(eqn, i - j)
end

function evaluateB(eqn::SpanwiseWH, z, i, j)
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

function evaluateF(eqn::SpanwiseWH, z, i)
    δ = wavenumber(eqn, 1)
    k₂ = wavenumber(eqn, 2)
    if isinf(z)
        return zero(z)
    else
        return -k₂ * interlacedcompliance(eqn, i - shift(eqn)) / (z + δ)
    end
end
