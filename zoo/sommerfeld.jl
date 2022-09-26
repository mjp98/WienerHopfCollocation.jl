
import WienerHopfCollocation: evaluateA, evaluateB, evaluateF

# ------------
# Wiener--Hopf equation
# ------------

struct SommerfeldWH{N,T,M} <: AbstractEquationWH{N,T}
    incident::PlaneWave{M,T}
    function SommerfeldWH(wave::PlaneWave{M,T}) where {M,T}
        return new{1,T,M}(wave)
    end
end

incident(eq::SommerfeldWH) = eq.incident
frequency(eq::SommerfeldWH) = eq |> incident |> frequency
wavevector(eq::SommerfeldWH) = eq |> incident |> wavevector

function evaluateA(problem::SommerfeldWH, z, i, j)
    k = frequency(problem)
    return sqrt(im * (z - k)) * sqrt(-im * (z + k))
end
function evaluateB(::SommerfeldWH, z, i, j)
    return one(z)
end
function evaluateF(problem::SommerfeldWH, z, i)
    δ, k₂, _ = wavevector(problem)
    return isinf(z) ? zero(z) : -k₂ / (z - δ)
end

# ------------
# Exact
# ------------

struct SommerfeldExact{N,T}
    incident::PlaneWave{N,T}
end

incident(sol::SommerfeldExact) = sol.incident
wavevector(sol::SommerfeldExact) = sol |> incident |> wavevector
frequency(sol::SommerfeldExact) = sol |> incident |> frequency

function (sol::SommerfeldExact{N,T})(α, u) where {N,T}
    kx, ky, kz = wavevector(sol)
    k = frequency(sol)
    γ = (α, u) -> u ? sqrt(-im * (α + k)) : sqrt(im * (α - k))
    if u
        return -γ(α, true) * (ky / (α - kx)) * (1 / γ(α, true) - 1 / γ(kx, true))
    else
        return -ky / (γ(kx, true) * (α - kx)) / γ(α, false)
    end
end
