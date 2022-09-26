##----------------
## Equation
##----------------


abstract type AbstractEquationWH{N,T} end

function evaluate_forcing(problem::AbstractEquationWH, z, i)
    return evaluateF(problem, z, i)
end

function evaluate_equation(problem::AbstractEquationWH, z, i, j)
    A = evaluateA(problem, z, i, j)
    B = evaluateB(problem, z, i, j)
    return A, B
end

function evaluate_forcing(problem::AbstractEquationWH, z)
    return evaluateF(problem, z)
end

function evaluate_equation(problem::AbstractEquationWH, z)
    A = evaluateA(problem, z)
    B = evaluateB(problem, z)
    return A, B
end

evaluateA(::AbstractEquationWH, z, i, j) = zero(z)
evaluateB(::AbstractEquationWH, z, i, j) = zero(z)
evaluateF(::AbstractEquationWH, z, i) = zero(z)

size(::AbstractEquationWH{N}, i) where {N} = N
size(eq::AbstractEquationWH) = (size(eq, 1), size(eq, 2))
eachindex(eq::AbstractEquationWH) = 1:size(eq, 1)
eachmatrixindex(eq::AbstractEquationWH) = Iterators.product(eachindex(eq), eachindex(eq))

function evaluateA(eq::AbstractEquationWH, z)
    return [evaluateA(eq, z, i, j) for (i,j) in eachmatrixindex(eq)]
end
function evaluateB(eq::AbstractEquationWH, z)
    return [evaluateB(eq, z, i, j) for (i,j) in eachmatrixindex(eq)]
end
function evaluateF(eq::AbstractEquationWH, z)
    return [evaluateF(eq, z, i) for i in eachindex(eq)]
end
