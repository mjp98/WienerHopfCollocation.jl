
using UnPack

struct HurdWH{N,T} <: AbstractEquationWH{N,T}
    k::Complex{T}
    μ::Complex{T}
    ν::Complex{T}
    θ₀::Complex{T}

    function HurdWH(k, μ, ν, θ₀)
        args = promote(complex(k), μ, ν, θ₀)
        new{2,real(eltype(args))}(args...)
    end
end

function evaluateA(eq::HurdWH, z, i, j)
    @unpack k, μ, ν = eq

    β = α -> sqrt(k - α) * sqrt(α + k)
    βz = β(z)

    (i == 1 && j == 1) && return -0.5 * ((βz + k * μ))
    (i == 1 && j == 2) && return -0.5 * (1 + k * μ * inv(βz))
    (i == 2 && j == 1) && return -0.5 * (-βz + k * ν)
    (i == 2 && j == 2) && return -0.5 * (1 + k * ν * inv(βz))
end

function evaluateB(::HurdWH, z, i, j)
    i == j && return one(z)
    return zero(z)
end

function evaluateF(eq::HurdWH, z, i)
    @unpack k, μ, ν, θ₀ = eq
    S₀ = sin(θ₀)

    i == 1 && return (k / (2im * π)) * (μ - S₀) / (z - k * cos(θ₀))
    i == 2 && return (k / (2im * π)) * (ν + S₀) / (z - k * cos(θ₀))
end
