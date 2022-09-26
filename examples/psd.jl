using WienerHopfCollocation
using StaticArrays
using Setfield

include("../zoo/planewave.jl")
include("../zoo/spanwise.jl")
include("../zoo/spanwise_gust.jl")

# Physical parameters
kvec = @SVector [3.346156743143974, 1.0, 0.0]
freq = 0.29266677053737383 + 0.0im
wave = PlaneWave{3,ComplexF64}(kvec, freq)

# Wiener Hopf problem
eqn = SpanwiseWH(wave, [1, -0.2, -0.2, 0.1, 0.1, 0.01, 0.02, 0.01, 0.02 + 0.1im], 9)


eqn = SpanwiseGustWH(freq, 20, kvec[3], [1, -0.2, -0.2, 0.1, 0.1, 0.01, 0.02, 0.01, 0.02 + 0.1im], 9)

prob = CommonSolve.init(ProblemWHC, eqn;
    space=Legendre(SqrtLine{-1 / 4}(0.0)),
    n=128, # number of basis functions
    m=256, # number of collocation points
    ifvanish=true, # enforce solution to vanish at ±∞
    ifsave=true, # load existing saved data, or compute and save
    solver=:qr    # Linear system factorization: lu or qr
)

# Solve collocation problem
data = precompute(prob);
sol = preallocate(prob);
Φ = CommonSolve.solve!(prob, sol, data)

# -----------
# Compute psd
# -----------

cache = cauchyfarfieldcache(prob, 100);

θ = π / 4

ff = LinRange(0.2, 2, 10)
M = 0.2


using ProgressMeter


@showprogress for f in ff
    prob = @set prob.equation.incident.ω = complex(float(f))
    prob = @set prob.equation.incident.k[1] = complex(float(f)) / M

    Φ = CommonSolve.solve!(prob, sol, data)

    Φ(f * cos(θ), true)[1]
end



function Π(k₁, k₃)
    T = 0.0025
    Lₜ = 1.2 * (2π)
    kₑ = sqrt(π) * SpecialFunctions.gamma(5 / 6) / (Lₜ * SpecialFunctions.gamma(1 / 3))
    k₁ₑ = k₁ / kₑ
    k₃ₑ = k₃ / kₑ
    scale = (4 * T^2) / (9π)
    k = (k₁ₑ^2 + k₃ₑ^2)
    ν = 7 / 3
    return scale * k / (1 + k)^ν
end

function powerspectraldensity(solution, cdata, fdata, M, θ)
    prob = solution.problem

    ret = 0.0

    N = (length(solution.ϕ)-1) ÷ 2
    for n in -N:N
        @set problem.equation.incident.wavevector[3] = n
        φ = solve!(problem, solution, cdata)
        k = frequency(problem)
        α =  -k * cos(θ)
        Ψ = φ(α , true, fdata)[n]
        δ = streamwise_wavenumber(problem)
        ret += abs2(Ψ) * Π(δ, n)

    end
    ϕ = ((1 - M * cos(θ)) * sin(θ))
    ret *= ϕ
    ret *= ϕ
    return ret
end

function computePSD(problem, ωs, M, θ)

    # Preallocate
    solution = preallocate(problem)
    # Precompute
    cdata = precompute(problem)
    # Use kwargs... so I remember what the arguments mean
    fdata = cauchyfarfieldcache(problem, 100, 1.0)

    PSD = Vector{Float64}(undef, length(ωs))

    @progress for (i, ω) in enumerate(ωs)

        @set problem.equation.incident.frequency = ω

        PSD[i] = powerspectraldensity(problem, solution, cdata, fdata, M, θ)

    end

    return PSD

end
