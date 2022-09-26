k = 1.0 + 0.0im
μ = sinpi(1 / 4) / (1im * k) #S1
ν = sinpi(1 / 5) / (1im * k) #S2
θ₀ = pi / 3

# Wiener Hopf problem
eqn = HurdWH(k, μ, ν, θ₀)
prob = CommonSolve.init(ProblemWHC, eqn)
Φ = CommonSolve.solve!(prob)
