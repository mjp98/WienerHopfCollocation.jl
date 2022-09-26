##--------------------------------------
## Solve wrapper
##--------------------------------------

function solve!(prob::ProblemWHC,
    sol::DataWHC=preallocate(prob),
    pre::CauchyData=precompute(prob)
)

    collocationsystem!(sol, prob, pre)

    if solver(prob) == :lu
        lusolve!(sol)
    elseif solver(prob) == :qr
        qrsolve!(sol)
    else
        @error "Implement $(solver(prob))"
    end

    return SolutionWH!(sol, prob)
end

function qrsolve!(sol::DataWHC{T}) where {T}
    if isbig(T)
        ldiv!(sol.x.blocks, qr!(sol.A.blocks), sol.b.blocks)
    else
        F = QR(decompose!(sol.ws, sol.A.blocks)...)
        ldiv!(sol.x.blocks, F, sol.b.blocks)
    end
end

function lusolve!(sol::DataWHC{T}) where {T}
    if isbig(T)
        ldiv!(sol.x.blocks, lu!(sol.A.blocks), sol.b.blocks)
    else
        F = LU(decompose!(sol.ws, sol.A.blocks)...)
        ldiv!(sol.x.blocks, F, sol.b.blocks)
    end
end
