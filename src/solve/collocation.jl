##--------------------------------------
## Collocation system
##--------------------------------------

function collocationsystem!(sol::DataWHC, prob::ProblemWHC, pre::CauchyData)
    collocationmatrix!(sol.A, prob, pre)
    collocationrhs!(sol.b, prob, pre)
end

function collocationmatrix(problem::ProblemWHC, cdata=cauchydata(problem))
    ret = preallocate_matrix(problem)
    return collocationmatrix!(ret, problem, cdata)

end

function collocationmatrix!(ret::AbstractMatrix, problem::ProblemWHC, cdata)
    if ifvanish(problem)
        return collocationmatrix_vanish!(ret, problem, cdata)
    else
        return collocationmatrix_novanish!(ret, problem, cdata)
    end
end

function collocationmatrix_novanish!(ret, problem::ProblemWHC, cdata)
    idx = 1:blocksize(problem, 1)

    for (I, J) in eachmatrixindex(problem)
        collocationblock!(ret, problem, cdata, I, J, idx)
    end

    return ret
end

function collocationmatrix_vanish!(ret, problem::ProblemWHC, cdata)
    idx = 2:blocksize(problem, 1)-1

    for (I, J) in eachmatrixindex(problem)

        if I == J
            vanish_firstrow!(ret, problem, I, J)
        end

        collocationblock!(ret, problem, cdata, I, J, idx)

        if I == J
            vanish_lastrow!(ret, problem, I, J)
        end

    end

    return ret
end

function collocationblock!(ret::AbstractMatrix, problem::ProblemWHC, cdata, I, J, idx)
    @unpack C₋, C₊, pts = cdata
    n = blocksize(problem, 2)

    for i in idx

        A, B = evaluate_equation(problem, pts[i], I, J)

        for j in 1:n

            ret[BlockIndex((I, J), (i, j))] = A * C₋[i, j] + B * C₊[i, j]
        end

    end
    return ret
end

function vanish_firstrow!(ret::AbstractMatrix, problem::ProblemWHC, I, J)
    n = blocksize(problem, 2)
    for j in 1:n
        ret[BlockIndex((I, J), (1, j))] = 1
    end
end

function vanish_lastrow!(ret::AbstractMatrix, problem::ProblemWHC, I, J)
    m, n = blocksize(problem)
    for j in 1:n
        ret[BlockIndex((I, J), (m, j))] = iseven(j) ? -1 : 1
    end
end


##--------------------------------------
## Collocation RHS
##--------------------------------------

function collocationrhs(problem::ProblemWHC)
    ret = problem |> preallocate_rhs
    x = problem |> collocationpoints
    return collocationrhs!(ret, problem, x)
end

function collocationrhs!(ret, problem::ProblemWHC, cdata::CauchyData)
    return collocationrhs!(ret, problem, collocationpoints(cdata))
end

function collocationrhs!(ret, problem::ProblemWHC, x::AbstractVector)
    if ifvanish(problem)
        return collocationrhs_vanish!(ret, problem, x)
    else
        return collocationrhs_novanish!(ret, problem, x)
    end
end

function collocationrhs_vanish!(ret, problem::ProblemWHC, x)
    n = length(x)
    for i in eachindex(problem)
        ret[BlockIndex((i), (1))] = 0
        for k in 2:(n-1)
            ret[BlockIndex((i), (k))] = evaluateF(problem, x[k], i)
        end
        ret[BlockIndex((i), (n))] = 0
    end
    return ret
end

function collocationrhs_novanish!(ret, problem::ProblemWHC, x)
    n = length(x)
    for i in eachindex(problem)
        for k in 1:n
            ret[BlockIndex((i), (k))] = evaluateF(problem, x[k], i)
        end
    end
    return ret
end
