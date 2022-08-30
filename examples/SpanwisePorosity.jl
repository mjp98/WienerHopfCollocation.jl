# add https://github.com/mjp98/HolomorphicFun.jl
# add https://github.com/mjp98/Aeroacoustics.jl
# add https://github.com/mjp98/RiemannHilbert.jl#update1_8
# add https://github.com/mjp98/WienerHopf.jl#update1_8
# add https://github.com/mjp98/ScaledLines.jl
# add https://github.com/mjp98/WienerHopfCollocation.jl


using Aeroacoustics
using CommonSolve
using HolomorphicFun
using WienerHopfCollocation
using Plots



numerical_options = Dict(
    :space => Legendre(SqrtLine{-1 / 4}(0.0)),
    # approximation space
    :n => 256, # number of basis functions
    :m => 256, # number of collocation points
    :ifvanish => true, # enforce solution to vanish at ±∞
    :ifsave => true, # load existing saved data, or compute and save
    :solver => :lu    # Linear system factorization: lu, qr, or sparse
);

incident = convectedgust(; M=0.1, ω=1)
variation = Fun(SawtoothWave(), 10)


todoubleinterlace(n::Int) = n >= 0 ? 2n + 1 : -2n
fromdoubleinterlace(n::Int) = iseven(n) ? -(n ÷ 2) : (n - 1) ÷ 2 # integer division

function getindexu(x::AbstractArray{T}, i::Int) where {T}
    j = todoubleinterlace(i)
    return j <= length(x) ? x[j] : zero(T)
end


function leadingedge_span(incident, compliance::Vector{T}) where {T}
    δ, k₂, k₃ = incident.k
    k = incident.ω
    N = length(compliance)
    w(n) = sqrt(k^2 - (k₃ + n)^2)
    c = compliance

    m = (N - 1) ÷ 2 # integer division

    # As elements
    A(z, i, j) = getindexu(c, i - j)
    function B(z, i, j)

        wₙ = sqrt(k^2 - (k₃ + j - m - 1)^2)

        g = sqrt(im * (z - wₙ)) * sqrt(-im * (z + wₙ))

        ret = getindexu(c, i - j) * g / 2

        if i == j
            ret += 1
        end
        return ret
    end
    F(z, i) = (-k₂ / (δ + z)) * getindexu(c, i - m - 1)

    # As matrix
    A(z) = T[A(z, i, j) for i = 1:N, j = 1:N]
    B(z) = T[B(z, i, j) for i = 1:N, j = 1:N]
    F(z) = T[F(z, i) for i = 1:N]

    return WienerHopfEq(A, B, F)
end


function LEspan_gust_problem(; incident, compliance, numerical_options)
    equation = leadingedge_span(incident, compliance)
    return WienerHopfCollocation.problem(equation, numerical_options)
end


spec = LEspan_gust_problem(; incident, compliance=coefficients(SawtoothWave(), 19), numerical_options)


cdata = precompute(spec);
ldata = preallocate(spec);
collocationmatrix!(ldata.A, spec, cdata);
collocationrhs!(ldata.b, spec, cdata);


# md"""Load packages"""
# end

# # ╔═╡ c7703f11-a103-4833-a9b6-88d2cad8eaff
# using RecursiveFactorization

# # ╔═╡ 46027a4f-f549-49aa-9a31-59040375504e
# using Setfield

# # ╔═╡ 71997923-9cff-44d7-b956-0464a922ec82
# using OffsetArrays

# # ╔═╡ 8e4fc0a7-3422-46a1-9ff5-ca1ab176ce4d
# using ProgressMeter

# # ╔═╡ bf67dc55-d950-4a5b-bbd5-138f4f3dc9e6
# md"""# Spanwise varying porosity"""

# # ╔═╡ 4901f998-677f-4f21-b698-f36bb947ec6f
# RHC = ingredients(path*"init_rewrite.jl");

# # ╔═╡ 22ccaf59-c6d5-45e1-83e5-d2ee0dce129c
# SqrtLineS{0.5}(complex(0.0),2.0)

# # ╔═╡ 06121bc4-b8be-4555-8712-ad2f678fe1e2
# md"""## Riemann--Hilbert collocation"""


# ╔═╡ cb22bcaf-a3e6-45f7-86df-7bde13c28e62
# begin
# 	struct AeroacousticDimensional{T}
# 		U::Quantity{T,SpeedDim}
# 		f::Quantity{T,FrequencyDim}
# 		k₃::Quantity{T,WavenumberDim}
# 	end

# 	struct AeroacousticNondimensional{T}
# 		M::T
# 		k::T
# 		k₃::T
# 	end

# 	function nondimensionalize(d::ReferenceDimensions,x::AeroacousticDimensional)
# 		return AeroacousticNondimensional(nondimensionalise(d,u) for u in x)
# 	end

# 	# In Prandtl--Glauert
# struct AeroacousticDisturbance{T}
# 	M::T
# 	k::Complex{T}
# 	wavevector::SVector{T}
# end

# 	struct PlateDimensional{T,S}
# 		ρ::Quantity{T,DensityDim}
# 		compliance::S
# 	end

# 	struct PlateNondimensional{T,S}
# 		ρ::T
# 		compliance::S
# 	end

# 	function nondimensionalize(d::ReferenceDimensions,x::PlateDimensional)
# 		return PlateNondimensional(nondimensionalise(d,x.ρ),x.compliance)
# 	end


# end

# ╔═╡ d65e526a-fd76-4e9a-bd93-4dfbeaaf5434
begin

    # 	function convectedgust(x::AeroacousticNondimensional)
    #  		β = sqrt(1-M^2)
    # 		δ = k/(β*M)
    #       wavevector = @SVector [δ,1,k₃]
    #       return AeroacousticDisturbance(M,k,wavevector)
    #   end

end

# ╔═╡ bc112860-8f14-4658-8f60-97cfb586bc9c
begin
    const WavenumberUnit = typeof(1.0u"1/m")
    const LengthUnit = typeof(1.0u"m")
    const SpeedUnit = typeof(1.0u"m/s")
    const FrequencyUnit = typeof(1.0u"1/s")
    const DensityUnit = typeof(1.0u"kg*m^-3")
    Base.@kwdef struct ThePhysics
        reference_length::LengthUnit = 0.1u"m"
        sound_speed::SpeedUnit = 343.0u"m/s"
        fluid_density::DensityUnit = 1u"kg*m^-3"
    end
    for element in [:reference_length, :sound_speed, :fluid_density]
        @eval $element(x::ThePhysics) = x.$element
    end
    reference_time(x::ThePhysics) = reference_length(x) / sound_speed(x)
    Base.length(x::ThePhysics) = reference_length(x)

    @with_kw mutable struct ParametersDim{T,S,W,X,V}
        physics::T
        α
        ρ::V = 1.0u"kg*m^-3"
        U::W = 50.0u"m/s"
        f::S = 500.0u"1/s"
        k₃::X = 0.0u"1/m"
    end


    nondim(c::ThePhysics, x::LengthUnit) = x / length(c)
    nondim(c::ThePhysics, x::WavenumberUnit) = x * length(c)
    nondim(c::ThePhysics, x::SpeedUnit) = x / sound_speed(c)
    nondim(c::ThePhysics, x::FrequencyUnit) = x * length(c) / sound_speed(c)
    nondim(c::ThePhysics, x::DensityUnit) = x / fluid_density(c)

    @with_kw struct ParametersNonDim{K,T,S,V}
        physics::K
        α
        ρ::V = 1.0
        M::S = 0.4
        k::T = one(ComplexF64)
        k₃::S = 0.0
    end

    function nondim(x::ParametersDim)
        ρ = nondim(x.physics, x.ρ)
        M = nondim(x.physics, x.U)
        k = nondim(x.physics, x.f)
        k₃ = nondim(x.physics, x.k₃)
        return ParametersNonDim(x.physics, x.α, ρ, M, k, k₃)
    end


    convectedgust(x::ParametersNonDim) = convectedgust(x.M, x.k, x.k₃)
    function convectedgust(M, k, k₃)
        β = sqrt(1 - M^2)
        δ = k / (β * M)
        wavevector = @SVector [δ, 1, k₃]
        return complex(δ * M), wavevector
    end

    @with_kw struct ParametersNonDim{K,T,S,V}
        physics::K
        α
        ρ::V = 1.0
        M::S = 0.4
        k::T = one(ComplexF64)
        k₃::S = 0.0
    end

    struct ParametersReduced{W,T,S,V}
        physics::W
        k::T
        wavevector::S
        compliance::V
    end
    function LeadingEdgeGust(x::ParametersNonDim)
        k, wavevector = convectedgust(x)
        compliance = x.ρ * x.α
        return ParametersReduced(x.physics, k, wavevector, compliance)
    end

    leadingedge_gust(x::ParametersDim) = LeadingEdgeGust(nondim(x))

end

# ╔═╡ 9fecaf25-b41a-43f4-b5f0-f14db0deabe4
begin
    const WavenumberDim = Unitful.𝐋^-1
    const LengthDim = Unitful.𝐋
    const SpeedDim = Unitful.𝐋 * Unitful.𝐓^-1
    const FrequencyDim = Unitful.𝐓^-1
    const DensityDim = Unitful.𝐌 * Unitful.𝐋^-3

    # NoUnits forces simplification of μm/m = 1e-6 and Hz*s = 1

    function nondimensionalize(d, x::Quantity{T,WavenumberDim,S}) where {T,S}
        return NoUnits(x * reference_length(d))
    end
    function nondimensionalize(d, x::Quantity{T,FrequencyDim,S}) where {T,S}
        return NoUnits(x * reference_time(d))
    end
    function nondimensionalize(d, x::Quantity{T,LengthDim,S}) where {T,S}
        return NoUnits(x / reference_length(d))
    end
    function nondimensionalize(d, x::Quantity{T,SpeedDim,S}) where {T,S}
        return NoUnits(x / reference_speed(d))
    end
    function nondimensionalize(d, x::Quantity{T,DensityDim,S}) where {T,S}
        return NoUnits(x / reference_density(d))
    end

    Base.@kwdef struct ReferenceDimensions{T}
        length::Quantity{T,LengthDim} = 0.1u"m"
        speed::Quantity{T,SpeedDim} = 343.0u"m/s"
        density::Quantity{T,DensityDim} = 1.0u"kg*m^-3"
    end

    reference_length(d::ReferenceDimensions) = d.length
    reference_speed(d::ReferenceDimensions) = d.speed
    reference_density(d::ReferenceDimensions) = d.density
    reference_time(d::ReferenceDimensions) = reference_length(d) / sound_speed(d)
end

# ╔═╡ 180d2704-ee9c-4e45-99c4-f5efcb794915
10.0u"cm" - 4u"mm"

# ╔═╡ e0f75939-d8ff-4d93-a658-9ee699ab6184
md"""###### Expansion in Fourier components"""

# ╔═╡ de71446b-eb1c-43b0-bf76-5fe71147240a
md"""
Compliance distribution: coefficients are interlaced, so ordered as

$(\ldots,x_{-1},x_0,x_1,\ldots) \to (x_0,x_{-1},x_1,\ldots)$
"""


# ╔═╡ 6d458921-e872-4a1b-b5ee-25ace6fc571b
md"""##### Wiener--Hopf equation"""

# ╔═╡ acd4c04c-aff5-4075-9849-9691f90cd505
md"""
$A(\alpha)\Psi_-(\alpha) + B(\alpha)\Psi_+(\alpha) = F(\alpha)$
"""

# ╔═╡ 25ec98fb-709a-4021-9e4e-a8b751638a95
# struct SpanwiseCompliantLeadingEdge{T} <: AbstractWienerHopfEquation
# 	A::Function
# 	B::Function
# 	F::Function
# 	param::T

# 	function SpanwiseCompliantLeadingEdge(reduced)

# 		k = reduced.k
# 		δ,k₂,k₃ = reduced.wavevector
# 		c = coefficients(reduced.compliance)

# 		m = (N-1)÷2 # integer division

# 		# Scalar entries
# 		A(z,i,j) = getindexu(c,i-j)

# 		function B(z,i,j)
# 			wₙ = sqrt(k^2 - (k₃ + j-m-1)^2)

# 			g = sqrt(im*(z-wₙ))*sqrt(-im*(z+wₙ))

# 			ret = getindexu(c,i-j)*g/2

# 			if i ==j
# 				ret+=1
# 			end
# 			return ret
# 		end

# 		F(z,i) = (-k₂/(δ+z))*getindexu(c,i-m-1)

# 		# Array valued
# 		A(z) = T[A(z,i,j) for i = 1:N,j=1:N]
# 		B(z) = T[B(z,i,j) for i = 1:N,j=1:N]
# 		F(z) = T[F(z,i) for i = 1:N]

# 		return
# 	end


# end




# ╔═╡ 88fd398f-a552-4b8b-b67f-7930627bce94
# begin
# 	abstract type AbstractWienerHopfEquation end

# maybe want indices(x) ...

# 	WienerHopfA(x,z) = T[WienerHopfA(x,z,i,j) for i = 1:truncation(x),j= 1:truncation(x)]
# 	WienerHopfB(x,z) = T[WienerHopfB(x,z,i,j) for i = 1:truncation(x),j= 1:truncation(x)]
# 	WienerHopfF(x,z) = T[WienerHopfF(x,z,i) for i = 1:truncation(x)]


# 	function WienerHopfEquation(x)
# 		A(z...) = WienerHopfA(x,z...)
# 		B(z...) = WienerHopfB(x,z...)
# 		F(z...) = WienerHopfF(x,z...)
# 		return WienerHopfEq(A,B,F)
# 	end

# 	Base.@kwdef struct SpanwiseCompliantLeadingEdge{T} <: AbstractWienerHopfEquation
# 		k::T
# 		wavevector
# 		compliance
# 		truncation
# 	end

# 	function SpanwiseCompliantLeadingEdge(p::ParametersReduced,N)
# 		return SpanwiseCompliantLeadingEdge(p.k,p.wavevector,p.compliance,2N-1)
# 	end

# 	ApproxFun.coefficients(x::SpanwiseCompliantLeadingEdge) = coefficients(x.compliance)
# 	wavevector(x::SpanwiseCompliantLeadingEdge) = x.wavevector
# 	truncation(x::SpanwiseCompliantLeadingEdge) = x.truncation

# 	wavenumber(x::SpanwiseCompliantLeadingEdge) = x.k
# 	streamwise_wavenumber(x::SpanwiseCompliantLeadingEdge) = x.wavevector[1]
# 	spanwise_wavenumber(x::SpanwiseCompliantLeadingEdge) = x.wavevector[3]

# 	maxmode(x::SpanwiseCompliantLeadingEdge) = (truncation(x)-1)÷2
# 	modes(x::SpanwiseCompliantLeadingEdge) = -maxmode(x):maxmode(x)

# 	shift(x::SpanwiseCompliantLeadingEdge) = (truncation(x)+1)÷2
# 	Base.size(x::SpanwiseCompliantLeadingEdge) = (truncation(x),truncation(x))

# 	function branchpoint(x::SpanwiseCompliantLeadingEdge,n)
# 		k  = x.k
# 		k₃ = x.wavevector[3]
# 		wₙ = sqrt(k^2 - (k₃ + n)^2)
# 		return wₙ
# 	end

# 	function scalarkernel(x::SpanwiseCompliantLeadingEdge,n::Int,z)
# 		wₙ = branchpoint(x,n)
# 		γₙ = sqrt(im*(z-wₙ))*sqrt(-im*(z+wₙ))
# 		return γₙ
# 	end

# 	function interlacedcoefficient(x::SpanwiseCompliantLeadingEdge,n::Int)
# 		return getindexu(coefficients(x),n)
# 	end

# 	function WienerHopfA(x::SpanwiseCompliantLeadingEdge,z,i,j)
# 		A = interlacedcoefficient(x,i-j)
# 		return A
# 	end

# 	function WienerHopfB(x::SpanwiseCompliantLeadingEdge,z,i,j)
# 		n  = j-shift(x)
# 		γₙ = scalarkernel(x,n,z)
# 		B = interlacedcoefficient(x,i-j)*γₙ/2
# 		if i == j
# 			B += 1
# 		end
# 		return B
# 	end

# 	function WienerHopfF(x::SpanwiseCompliantLeadingEdge,z,i)
# 		δ,k₂,_ = x.wavevector
# 		F = -k₂*interlacedcoefficient(x,i-shift(x))/(δ+z)
# 		return F
# 	end

# 	function velocitypotential(Ψ::SolutionWHC{SpanwiseCompliantLeadingEdge},α)
# 		_,_,k₃ = wavevector(Ψ)
# 		Ψ₊ = Ψ(α,true)
# 		ψ = complex(0.0)
# 		for m = modes(Ψ)
# 			ret += Ψ₊[m]*exp(im*(k₃+m)*z - scalarkernel(Ψ,m,α)*abs(y))
#		end
#		return sign(y)*ret/2
# 		end
# 		return Φ
# 	end

# struct Point3DCyl{T<:Real}
# 	r::T
# 	θ::T
#   z::T
# end

# struct Point2DCyl{T<:Real}
# 	r
# 	θ
# end

# 	function farfieldpressuremode(Ψ,n::Int,x::Point3DCyl)
# 		r,θ,z = x
# 		δ,_,k₃ = wavevector(Ψ)
# 		wₙ = branchpoint(Ψ,n)
# 		α = -wₙ*cos(θ)
# 		Ψₙ₊ = Ψ(α,true)[n]
# 		scale = cispi(1/4)*sin(θ)/sqrt(2π*r)
# 		return scale*Ψₙ₊*sqrt(wₙ)*exp(im*(wₙ*r + (k₃ + n)*z))*(δ+α)
# 	end

# 	modes(x::SolutionWHC{SpanwiseCompliantLeadingEdge}) = modes(x.equation)

# 	function farfieldpressure(Ψ::SolutionWHC{SpanwiseCompliantLeadingEdge},x::Point3DCyl)
# 		p = complex(0.0)
# 		for m ∈ modes(Ψ)
# 			p += farfieldpressuremode(Ψ,n,x)
# 		end
# 		return p
# 	end

# OffsetArray(Ψ,-N:N)

# end

# ╔═╡ 6923f795-8400-42bc-8fcb-6eff1e48fbc3
begin

    struct Waveform <: Function
        f
        coefficients
        name::Symbol
    end
    (f::Waveform)(x) = f.f(x)
    ApproxFun.coefficients(f::Waveform, n) = f.coefficients(n)
    ApproxFun.Fun(f::Waveform, N::Integer) = Fun(Fourier(PeriodicSegment(0 .. 1)), [coefficients(f, n) for n = 0:N])

    function sawtooth(t)
        return t - floor(t)
    end
    function sawtoothc(n)
        if n == 0
            return 1 / 2
        elseif mod(n, 2) == 1
            m = div(n + 1, 2)
            return -1 / (m * pi)
        else
            return 0
        end
    end

    SawtoothWave = Waveform(sawtooth, sawtoothc, :sawtooth)

end

# ╔═╡ 6043ac69-567c-4544-8cf7-092f232e4a20
begin
    periodicvariation(f::Function) = Fun(z -> f(z), Laurent(PeriodicSegment(0 .. 2π)))
    periodicvariation(f::Function, m) = Fun(z -> f(z), Laurent(PeriodicSegment(0 .. 2π)), m)
end

# ╔═╡ 907d8928-be18-4e56-8c06-3e8a0e653277
begin
    leadingedge_span(p::ParametersReduced, n) = leadingedge_span(p.k, p.wavevector, p.compliance, n)

    function leadingedge_span(k::T, wavevector, compliance, N::Int) where {T}
        δ, k₂, k₃ = wavevector
        w(n) = sqrt(k^2 - (k₃ + n)^2)

        c = coefficients(compliance)

        m = (N - 1) ÷ 2 # integer division

        # As elements
        A(z, i, j) = getindexu(c, i - j)
        function B(z, i, j)

            wₙ = sqrt(k^2 - (k₃ + j - m - 1)^2)

            g = sqrt(im * (z - wₙ)) * sqrt(-im * (z + wₙ))

            ret = getindexu(c, i - j) * g / 2

            if i == j
                ret += 1
            end
            return ret
        end
        F(z, i) = (-k₂ / (δ + z)) * getindexu(c, i - m - 1)

        # As matrix
        A(z) = T[A(z, i, j) for i = 1:N, j = 1:N]
        B(z) = T[B(z, i, j) for i = 1:N, j = 1:N]
        F(z) = T[F(z, i) for i = 1:N]

        return WienerHopfEq(A, B, F)
    end

end

# ╔═╡ 7951dc15-62e5-4497-a20c-c4885726f69b
function LEspan_gust_problem(dimensional, truncation, numerical_options)
    nondim_gust = leadingedge_gust(dimensional)
    equation = leadingedge_span(nondim_gust, truncation)
    return problem(equation, numerical_options)
end

function LEspan_gust_problem(dimensional, truncation, numerical_options)
    equation = leadingedge_span(nondim_gust, truncation)
    return problem(equation, numerical_options)
end

# ╔═╡ 6a63b65f-b289-439f-90f7-8b5f561b241a
function baselineparam(x::ParametersDim)
    @unpack physics, α, ρ, U, f, k₃, = x
    return ParametersDim(physics, Fun(space(α), [coefficients(α)[1]]), ρ, U, f, k₃)
end

# ╔═╡ dfe8e412-2269-432f-bc6f-8370205d49db
function baselineparam(x::ParametersReduced)
    meancompliance = Fun(space(x.compliance), [coefficients(x.compliance)[1]])
    return newx = @set x.compliance = meancompliance
end

# ╔═╡ 0d65e2cb-204f-4eab-bcc7-8e0d693013d1
F = Fun(SawtoothWave, 10)

# ╔═╡ ec8809d6-6347-4909-89a8-51e88bcb5f5e
G = Fun(z -> 1 + (F(z) - 0.5), Laurent(PeriodicSegment(0 .. 1)))

# ╔═╡ 9895911a-aaf2-4f20-81cf-bbdbebc6c20a
begin

    numerical_options = Dict(
        :space => Legendre(ScaledLine(SqrtLine{-1 / 4}(0.0), 1)),
        # approximation space
        :n => 256, # number of basis functions
        :m => 256, # number of collocation points
        :ifvanish => true, # enforce solution to vanish at ±∞
        :ifsave => true, # load existing saved data, or compute and save
        :solver => :lu    # Linear system factorization: lu, qr, or sparse
    )

    physical_constants = ThePhysics(
        sound_speed=343.0u"m/s",
        reference_length=36u"mm"
    )

    dimensional = ParametersDim(
        physics=physical_constants,
        α=periodicvariation(z -> (1 + 0.9 * cos(z))),
        f=500.0u"1/s",
        U=0.2 * 343.0u"m/s",
        k₃=0.00u"1/m"
    )

    N = 5
end;

# ╔═╡ a8a8a7dc-c2ae-4ce9-97c1-a20ae3a34620
begin
    plotcompliance(x::ParametersDim; kwargs...) = plotcompliance(x.α; kwargs...)
    plotcompliance(x::ParametersReduced; kwargs...) = plotcompliance(x.compliance; kwargs...)
    function plotcompliance(C::Fun; kwargs...)
        plot(C; lw=2, labels=["real" "imag"])
        plot!(xlims=(0, 2π), legend=:outertopright, ylabel=L"C(z)", xlabel=L"z")
        plot!(; kwargs...)
    end
end

# ╔═╡ 60b92a45-7ead-45fe-9f1c-fcbf6d6df003


# ╔═╡ f414a231-4dcb-44f5-9af7-5084d6c29594


# ╔═╡ 60292092-d294-4bd1-aea4-d8f64ac875fb
plotcompliance(dimensional; title="compliance variation")

# ╔═╡ df9278ab-7072-41c1-95d9-788e537ed311
md"""### Uniform porosity"""

# ╔═╡ 513eddf0-9475-4b48-b655-d113f11a6e5e
# begin
# function evaluate(solution::SolutionWHC{SpanwiseCompliantEdge},z::T) where T<:Number
#     inUHP = isabove(space(solution),z)
#     cz = cauchy(solution,z)
#     @unpack A,B,F = equation(solution)
#     if inUHP
#         u₊ = cz
# 		u₋ = A(z) \ (F(z) - B(z) * cz)
#     else
#         u₊ = B(z) \ (F(z) - A(z) * cz)
# 		u₋ = cz
#     end
# 	idx = modes(solution)
# 	return OffsetArray(u₊,idx),OffsetArray(u₋,idx)
# end


# function evaluate(solution::SolutionWHC{SpanwiseCompliantEdge},z::T,isUHP::Bool) where T<:Number
#     @unpack A, B, F = equation(solution.problem)
#     inUHP = isabove(solution,z)
#     cz = cauchy(solution,z)
#     if isUHP
#         u = inUHP ? cz : B(z) \ (F(z) .- A(z) * cz)
#     else
#         u = inUHP ? A(z) \ (F(z) .- B(z) * cz) : cz
#     end
# 	idx = modes(solution)
# 	return OffsetArray(u,idx)
# end
# end

# ╔═╡ 8a459409-d622-47dc-947b-b8ee82791c4f
begin
    struct PhysicalSolution{T<:ParametersReduced,S<:AbstractSolutionWHC}
        parameters::T
        φ::S
    end
    Base.size(x::PhysicalSolution) = size(x.φ)
    (x::PhysicalSolution)(z...) = x.φ(z...)
    evaluatecached(z, x...) = evaluatecached(z.φ, x...)

    function velocitypotential(Ψ::PhysicalSolution)
        (N, _) = size(Ψ)
        k = Ψ.parameters.k
        wavevector = Ψ.parameters.wavevector
        _, _, k₃ = wavevector
        M = (N - 1) ÷ 2
        @inline w(n) = sqrt(k^2 - (k₃ + n)^2)
        function Φ(α, y, z)
            ret = zero(α)
            Ψ₊ = Ψ(α, true)
            for m = -M:M
                ret += Ψ₊[m+M+1] * exp(im * (k₃ + m) * z - helmholtzγ(α, w(m)) * abs(y))
            end
            return sign(y) * ret / 2
        end
        return Φ
    end

    function farfieldpressure(Ψ::PhysicalSolution)
        (N, _) = size(Ψ)
        k = Ψ.parameters.k
        wavevector = Ψ.parameters.wavevector
        δ, k₂, k₃ = wavevector
        M = (N - 1) ÷ 2
        function p(r, θ::T, z) where {T}
            ret = zero(complex(T))
            for m = -M:M
                wm = sqrt(k^2 - (k₃ + m)^2)
                α = -wm * cos(θ)
                Ψ₊ = Ψ(α, true)[m+M+1]
                ret += Ψ₊ * sqrt(wm) * exp(im * (wm * r + (k₃ + m)) * z) * (δ + α)
                #	end
            end
            return cispi(1 / 4) * sin(θ) * ret / sqrt(2π * r)
        end
        return p
    end

    function farfieldpressure(Ψ::PhysicalSolution, cache)
        (N, _) = size(Ψ)
        k = Ψ.parameters.k
        wavevector = Ψ.parameters.wavevector
        δ, k₂, k₃ = wavevector
        M = (N - 1) ÷ 2
        function p(r, θ::T, z) where {T}
            ret = zero(complex(T))
            for m = -M:M
                wm = sqrt(k^2 - (k₃ + m)^2)
                α = -wm * cos(θ)
                Ψ₊ = evaluatecached(Ψ, α, true, cache)[m+M+1]

                # If we integrate over z, then we just need integral 0 to 2π of exp(im*(k₃ + m)*z) ???

                ret += Ψ₊ * sqrt(wm) * exp(im * (wm * r + (k₃ + m) * z)) * (δ + α)
            end
            return cispi(1 / 4) * sin(θ) * ret / sqrt(2π * r)
        end
        return p
    end

end

# ╔═╡ 2bc7e143-8910-453f-816c-afd8ea30934c
begin
    spec = LEspan_gust_problem(dimensional, N, numerical_options)

    precomputed = precompute(spec)# collocation points and Cauchy matrices
    preallocated = preallocate(spec)# linear system and solution vector

    φ = lusolve!(preallocated, precomputed, spec)

    ϕ = PhysicalSolution(leadingedge_gust(dimensional), φ)
end;

# ╔═╡ 52088b36-b3e3-4fef-b1a1-3c59da0a3261
plotequation(φ; title="block sparsity")

# ╔═╡ 6ebe0310-e3a0-47ce-885c-37fbc545637a
plotcoefficients(φ; cs=:turbo, ylims=(1e-20, 100), xlims=(0, 200), legend=:outertopright, grid=:none)

# ╔═╡ 1608076e-0503-440f-a993-8e32b66dca92
begin
    function branchcut_endpoints(wₙ; ymax=10)
        yₙ = max(ymax, abs(imag(wₙ)))
        x = real(wₙ) * [1, 1]
        if imag(wₙ) < 0 || (imag(wₙ) == 0 && real(wₙ) < 0)
            yₙ *= -1
        end
        y = [imag(wₙ), yₙ]
        return x, y
    end

    function plot_branchcuts(k, k₃, m; xscale=3, yscale=3, kwargs...)
        N = (m - 1) ÷ 2
        cs = [get(ColorSchemes.colorschemes[:turbo], i) for i in LinRange(0, 1, 2 * N + 1)]
        plt = plot(palette=cs)
        ymax = yscale * abs(k)
        xmax = xscale * abs(k)
        for n = -N:N
            wₙ = sqrt(complex(k^2 - (k₃ + n)^2))

            x, y = branchcut_endpoints(wₙ; ymax)
            scatter!([x[1]], [y[1]], c=cs[n+N+1], markersize=2, label=:none)
            scatter!(-[x[1]], -[y[1]], c=cs[n+N+1], markersize=2, label=:none)
            plot!(-x, -y; label=:none, c=cs[n+N+1], kwargs...)
            plot!(x, y; label=latexstring(n), c=cs[n+N+1], kwargs...)
        end
        plot!(ylims=(-ymax, ymax), xlims=(-xmax, xmax), legend=:outertopright, grid=:none, framestyle=:origin)
        return plt
    end

    function plot_branchcuts(ϕ::PhysicalSolution; kwargs...)
        k = ϕ.parameters.k
        k₃ = ϕ.parameters.wavevector[3]
        m, _ = size(ϕ.φ.problem.equation)
        plot_branchcuts(k, k₃, m; kwargs...)
    end
end


# ╔═╡ beaac2b0-f1d2-45df-b0aa-4cad4404cdcc
# if k₃ = 0 then negative and positive branch cuts coincide
plot_branchcuts(ϕ; lw=2, xscale=2, yscale=100, title="branch cuts")

# ╔═╡ e741c833-8d69-44c8-933b-534d9a04a72e
begin
    baseparam = baselineparam(dimensional)
    basespec = LEspan_gust_problem(baseparam, 1, numerical_options)
    φbase = solve(basespec)
    ϕbase = PhysicalSolution(leadingedge_gust(baseparam), φbase)
end;

# ╔═╡ f34c07b1-3e5a-46d3-add9-392dcb3876da
plotcompliance(baseparam)

# ╔═╡ 84d1a531-8db2-4563-998a-f27c4af1954b
fcache = cauchyfarfieldcache(spec, 400, 10.0);

# ╔═╡ 5485c4bb-b78b-4cc1-9411-870045caacca
begin
    θd = LinRange(0.0, 2π, 300)
    r = 10
end;

# ╔═╡ fc907d25-95ff-4cdc-827c-07a5f9b678f8
md"""### Tuncation convergence"""

# ╔═╡ 1b3cbd7c-6915-4b7a-aafc-ee2e4c98d2c4
function smallerdata!(sdata, ldata)
    (M, _) = blocksize(ldata.A)
    (N, _) = blocksize(sdata.A)

    Δ = (M - N) ÷ 2

    for i = 1:N
        sdata.x[Block(i)] .= ldata.x[Block(i + Δ)]
        sdata.b[Block(i)] .= ldata.b[Block(i + Δ)]
        for j = 1:N
            sdata.A[Block(i, j)] .= ldata.A[Block(i + Δ, j + Δ)]
        end
    end
    return sdata
end

# ╔═╡ 5faab6ae-1455-4355-9936-403ee67bd86b
maxtrunc = 21;

# ╔═╡ fb7f308e-2f41-4134-9a0d-df383660b5b5
begin
    maxspec = LEspan_gust_problem(dimensional, maxtrunc, numerical_options)
    cdata = precompute(maxspec)
    ldata = preallocate(maxspec)
    collocationmatrix!(ldata.A, maxspec, cdata)
    collocationrhs!(ldata.b, maxspec, cdata)

    trad = 2.0
    tθ = π / 3
end;

# ╔═╡ 0b36d1a2-fc8c-48c0-bcf5-f533bfd43c29
function truncationplot!(plt, ps)
    for (i, rs) in enumerate(ps[1:end-1])
        plt = scatter!([2i - 1], [abs.(1 .- rs / ps[end])]; label=:none, c=:black, markersize=2)
    end
    return plt
end

# ╔═╡ 6ebd839d-acef-4b97-ac9d-d77fcc13ce09
begin
    reduced = leadingedge_gust(dimensional)
    @set reduced.wavevector[3] = 2.0
end

# ╔═╡ 1bac3865-a51b-4b6b-8738-467404290de5
function Π(k₁, k₃)
    T = 0.0025
    Lₜ = 1.2 * (2pi)
    kₑ = sqrt(pi) * SpecialFunctions.gamma(5 / 6) / (Lₜ * SpecialFunctions.gamma(1 / 3))
    k₁ₑ = k₁ / kₑ
    k₃ₑ = k₃ / kₑ
    scale = (4 * T^2) / (9π)
    k = (k₁ₑ^2 + k₃ₑ^2)
    ν = 7 / 3
    return scale * k / (1 + k)^ν
end

# ╔═╡ 9d7ff067-e79e-47e7-83bf-00b94fd93bf3

function computePSD(freqrange, dim, ex, r, θ, z)

    # setup
    preall = preallocate(ex)
    precom = precompute(ex)
    fcache = cauchyfarfieldcache(ex, 100, 1.0)

    M = nondim(dim).M

    PSD = Vector{Float64}(undef, length(freqrange))

    nwh, _ = size(ex)
    N = (nwh - 1) ÷ 2
    opts = numerics(ex)

    #println("Computing PSD")
    @progress for (i, f) in enumerate(freqrange)
        # update frequency
        newdim = @set dim.f = f
        # reduce parameters
        reduced = LeadingEdgeGust(nondim(newdim))
        # compute psd
        psd = 0.0
        for n = -N:N
            @set reduced.wavevector[3] = n
            #problem = LEspan_gust_problem(reduced,nwh,opts)

            equation = leadingedge_span(reduced, nwh)
            problem = ProblemWHC(equation, opts)

            φ = solve!(preall, precom, problem)
            Ψₙ₊ = evaluatecached(φ, -reduced.k * cos(θ), true, fcache)[n+N+1]
            #Ψₙ₊=1
            δ = reduced.wavevector[1]

            psd += abs2(Ψₙ₊) * Π(δ, n)
        end
        ϕ = ((1 - M * cos(θ)) * sin(θ))
        psd *= ϕ
        psd *= ϕ

        PSD[i] = psd
    end
    return PSD
end


# ╔═╡ cb02c6b7-79bf-4b61-9532-56b2979a19c3
freqrange = (10 .^ LinRange(log10(20), log10(20000), 140)) * 1.0u"1/s"

# ╔═╡ b51b631d-e233-4382-93a7-e89f90109238


# ╔═╡ d60b3d2d-f2ac-44dc-9326-98449ba5efd0
psd = computePSD(freqrange, dimensional, deepcopy(spec), 1.0, π / 4, 0.0)

# ╔═╡ c4604a10-f74f-4060-ab62-29650f4cc38c
psdbase = computePSD(freqrange, baseparam, deepcopy(spec), 1.0, π / 4, 0.0)

# ╔═╡ 907b5082-d752-4af6-950a-b0209a8f004b
begin
    plot()
    plot!(ustrip.(freqrange), 10 * log10.(psdbase); label="uniform", lw=2, xscale=:log10)
    plot!(ustrip.(freqrange), 10 * log10.(psd); label="vary", lw=2, xscale=:log10, ls=:dash)
    plot!(xlabel="frequency [Hz]", ylabel="PSD [dB]")
end

# ╔═╡ 59d521ac-719d-4ac5-96e1-37b5db01ca5e
begin
    plot()
    plot!(ustrip.(freqrange), 10 * log10.(psd ./ psdbase); label=:none, lw=2, xscale=:log10, ylims=(-2, 2), c=:black)
    plot!(xlabel="frequency [Hz]", ylabel=" change PSD [dB]")
end

# ╔═╡ 775b94cb-4a6a-4010-ba04-06e6a173d2d1
begin
    physics = ThePhysics()

end

# ╔═╡ d71a3402-5c9c-4213-a2ae-0cac3c0cb741
nondimensionalize(physics, 0.4u"cm")

# ╔═╡ dd728841-e069-402a-9e0e-ed4869105f5d
function reduced_parameters(; physics=ThePhysics(),
    M=0.4,
    k₁=2,
    k₂=1.0,
    k₃=0.0,
    massdistribution=z -> 0.1
)

    compliance = Fun(massdistribution, Laurent(PeriodicSegment(0, 2π)))
    β = sqrt(1 - M^2)
    δ = k₁ / β
    k = complex(δ * M)
    wavevector = @SVector [δ, 1.0, k₃]
    return ParametersReduced(physics, k, wavevector, compliance)
end


# ╔═╡ df910db8-0a30-43e9-8690-69c1f6bb2957
begin
    a = 0.45
    m₀ = 1
    m₀ = 0.698
    #m₀ = 0.327
    massdistribution = z -> m₀ * (1 + 2a * cos(z))

    reducedAIAA = reduced_parameters(; M=0.4, k₁=0.11, massdistribution=massdistribution)
end

# ╔═╡ b5acbfd9-1bbc-40c9-b955-4084d8bdd020
function farfieldpressure(x::ParametersReduced, truncation, opts, fcache)
    equationAIAA = leadingedge_span(reducedAIAA, truncation)

    problemAIAA = problem(equationAIAA, opts)

    φAIAA = solve(problemAIAA)
    ϕAIAA = PhysicalSolution(reducedAIAA, φAIAA)

    pAIAA = farfieldpressure(ϕAIAA, fcache)
    return pAIAA
end

# ╔═╡ c9a14010-db1b-4d3d-ac7a-05f960466ae4
pt = farfieldpressure(ϕ, fcache);

# ╔═╡ 445c5108-f4fe-43ab-8fd7-009bfd5bcb75
PlutoUI.with_terminal() do
    @code_warntype pt(10.0, pi / 3, 0.0)
end

# ╔═╡ 58bdafaf-86cf-4e0c-912b-01bd7077be51
pb = farfieldpressure(ϕbase, fcache);

# ╔═╡ 721b33ff-1505-4654-b277-2762e1ac1ac0
begin
    ptest = [pt(r, θ, 0.0) for θ in θd]
    pbase = [pb(r, θ, 0.0) for θ in θd]
end

# ╔═╡ 0d5c6fc6-ef13-4e19-8391-002986cb461b
begin
    plot()
    plot!(θd, abs.(pbase); proj=:polar, lw=2, label="uniform")
    plot!(θd, abs.(ptest); proj=:polar, lw=2, label="vary", ls=:dash)
    plot!(legend=:outertopright)
end

# ╔═╡ 49ba2cb2-5642-4189-8610-282bf40c35f7
begin
    nsmaller = 1:2:maxtrunc
    ps = Vector{ComplexF64}(undef, length(nsmaller))
    tparam = leadingedge_gust(dimensional)

    @progress for (i, t) in enumerate(nsmaller)

        tspec = LEspan_gust_problem(dimensional, t, numerical_options)
        tdata = preallocate(tspec)

        smallerdata!(tdata, ldata)
        sparsesolve!(tdata)

        tφ = SolutionWHC!(tdata.φ, tspec, tdata.x)
        tϕ = PhysicalSolution(tparam, tφ)

        tp = farfieldpressure(tϕ, fcache)
        ps[i] = tp(trad, tθ, 0)
    end
end

# ╔═╡ 36a2db17-5bdf-4531-ada3-4e211fb61390
begin
    plt = plot(grid=:none)
    plt = truncationplot!(plt, ps)
    plt = plot!(title="far-field pressure: truncation convergence")
    plt = plot!(xlabel=L"N", ylabel="relative error")
    plt = plot!(ylims=(1e-15, 1), yscale=:log10)
end

# ╔═╡ e01739d3-62e5-4872-8f7d-65d749b90580
reducedAIAA0 = baselineparam(reducedAIAA)

# ╔═╡ 1da2d1bd-e99a-4aa7-9bef-0e13be8eb433
plotcompliance(reducedAIAA0)

# ╔═╡ c397a08d-23bb-4da2-b69a-d2e0e531d542
plotcompliance(reducedAIAA)

# ╔═╡ 564dae62-00f9-4fd2-a4a7-97348245b490
pAIAA = farfieldpressure(reducedAIAA, 11, numerical_options, fcache);

# ╔═╡ 9388dba6-016f-40be-a23b-8c89ee781159
pAIAA0 = farfieldpressure(reducedAIAA0, 1, numerical_options, fcache);

# ╔═╡ 7873cbeb-5a8c-4507-ba92-411c8f90bf3d
begin
    pzAIAA = [pAIAA(r, θ, 0.0) for θ in θd]
    pzAIAA0 = [pAIAA0(r, θ, 0.0) for θ in θd]
end

# ╔═╡ c3625e88-f862-4854-9f62-55c8e2e44cab
begin
    plot()
    plot!(θd, abs.(pzAIAA0); proj=:polar, lw=2, label="uniform")
    plot!(θd, abs.(pzAIAA); proj=:polar, lw=2, label="cos", ls=:dash)
end

# ╔═╡ 7862d882-d7ae-45a1-9b66-7cc58d02625d
plotcompliance(reducedAIAA0)

# ╔═╡ 52504e10-6602-49ac-bb01-13d11e178561
plotcompliance(reducedAIAA)
