module WienerHopfCollocation

using Reexport

@reexport using ApproxFun
@reexport using CommonSolve
@reexport using ScaledLines
@reexport using SingularIntegralEquations
@reexport using WienerHopf

using BlockArrays
using FastLapackInterface
using JLD2
using Lazy
using LinearAlgebra
using OffsetArrays
using RiemannHilbert
using StaticArrays
using UnPack

import ApproxFun: prectype, coefficients
import Base: size, eachindex, eltype, length, sign
import BlockArrays: blocksizes, blockaxes
import CommonSolve: init, solve!
import RiemannHilbert: collocationpoints, fpstieltjesmatrix
import RiemannHilbert: evaluationmatrix!, evaluationmatrix
import RiemannHilbert: orientedleftendpoint, orientedrightendpoint
import SingularIntegralEquations: cauchy, stieltjesmoment!

export isabove, preallocate, precompute
export AbstractEquationWH, ProblemWHC, OptionsWHC, SolutionWH!
export collocationrhs, collocationrhs!
export collocationmatrix, collocationmatrix!
export qrsolve!, lusolve!
export cauchyfarfieldcache, CauchyCache
export evaluate_equation, evaluate_forcing


RiemannHilbert.orientedleftendpoint(d::ScaledSegmentLine) = RiemannDual(leftendpoint(d.line), RiemannHilbert.sign(d.line))
RiemannHilbert.orientedrightendpoint(d::ScaledSegmentLine) = RiemannDual(rightendpoint(d.line), RiemannHilbert.sign(d.line))

include("util/util.jl")
include("structure/structure.jl")
include("solve/solve.jl")

end
