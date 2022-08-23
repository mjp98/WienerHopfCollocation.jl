##--------------------------------------
## Split function helper
##--------------------------------------

abstract type AbstractSplit <: Function end

abstract type AbstractSumSplit <: AbstractSplit end
abstract type AbstractProdSplit <: AbstractSplit end
abstract type AbstractSolutionSplit <: AbstractSplit end

abstract type AbstractFactor end

struct ProdFactor{T,S} <: AbstractFactor
    x::S
end

struct SumFactor{T,S} <: AbstractFactor
    x::S
end

struct SolutionFactor{T,S} <: AbstractFactor
    x::S
    function SolutionFactor{T}(x::S) where {T,S}
        return new{T,S}(x)
    end
end

(x::ProdFactor{T})(z) where T = evaluate(x.x,z,T)
(x::SumFactor{T})(z) where T = evaluate(x.x,z,T)
(x::SolutionFactor{T})(z) where T = evaluate(x.x,z,T)
