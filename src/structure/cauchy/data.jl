struct CauchyData{T}
    C₋::Matrix{T}
    C₊::Matrix{T}
    pts::Vector{T}
end

collocationpoints(x::CauchyData) = x.pts

precompute(problem) = cauchydata(problem)

function cauchydata(problem)
    if ifsave(problem)
        C₋, C₊, pts = cachedcauchydata(problem)
    else
        C₋, C₊, pts = newcauchydata(problem)
    end
    return CauchyData(C₋, C₊, pts)
end

function newcauchydata(problem)
    pts = collocationpoints(problem)
    C₋, C₊ = cauchymatrices(problem, pts)
    return C₋, C₊, pts
end

function cachedcauchydata(problem)
    filename = cauchydata_filename(problem)
    if isfile(filename)
        C₋, C₊, pts = load(filename, "cauchy_minus", "cauchy_plus", "pts")
    else
        C₋, C₊, pts = newcauchydata(problem)
        save(filename, "cauchy_minus", C₋, "cauchy_plus", C₊, "pts", pts)
    end
    return C₋, C₊, pts
end

function cauchydata_filename(problem)
    sp = space(problem)
    m = ncollocation(problem)
    n = nbasis(problem)
    return cauchydata_filename(sp, m, n)
end

function cauchydata_filename(sp, m, n)
    !isdir("cauchydata") && mkdir("cauchydata")

    filename = string("cauchydata/m$(m)_n$(n)")

    T = prectype(sp)

    if domain(sp) isa ScaledSegmentLine
        filename *= string("scale_$(ScaledLines.scale(sp))")
    end

    if isbig(T)
        filename *= string("_oftype_$(T)_$(precision(T)).jld2")
    else
        filename *= string("_oftype_$(T).jld2")
    end

    return filename
end







function cauchymatrices(problem::ProblemWHC, pts)
    C₋ = fpcauchymatrix(problem)
    E  = evaluationmatrix(problem, pts)
    C₊ = C₋ + E
    return C₋, C₊, pts
end

collocationpoints(x::ProblemWHC) = collocationpoints(space(x), collocation(x))
collocationpoints(x::Space, s::CollocationSpec) = collocationpoints(x, ncollocation(s))

function fpstieltjesmatrix(problem::ProblemWHC)
    sp = space(problem)
    m = ncollocation(problem)
    n = nbasis(problem)
    return fpstieltjesmatrix(sp, domain(sp), m, n)
end

function evaluationmatrix(problem::ProblemWHC, pts)
    return evaluationmatrix(basis(problem), pts)
end

function evaluationmatrix(basis::BasisSpec, pts::Vector{T}) where {T}
    @unpack n, sp = basis
    m = length(pts)
    E = Matrix{T}(undef, m, n)
    return evaluationmatrix!(E, sp, copy(pts))
end
