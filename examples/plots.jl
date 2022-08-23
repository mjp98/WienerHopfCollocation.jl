
plotequation(solution::SolutionWHC;kwargs...) = plotequation(solution.problem;kwargs...)
plotequation(problem::ProblemWHC;kwargs...) = plotequation(problem.equation;kwargs...)

function plotequation(equation::WienerHopfEq;title="Linear system sparsity",kwargs...)
	heatmap(log10.(abs.([Array(equation.A(0))+Array(equation.B(0)) Array(equation.F(0))]));yflip=true,c=:Greys,colorbar=:none,axis=:off,ticks=:none,aspectratio=1,title)
end

function plotcoefficients(Φ;cs=:turbo,yscale=:log10,ylims=(1e-15,1),xlims=(0,250),kwargs...)
    n = length(Φ.ϕ)
 	plot(palette=[get(ColorSchemes.colorschemes[cs],i) for i = LinRange(0,1,n)])
    index =  -Int((n-1)/2):Int((n-1)/2)
	m = (n-1)÷2

	labels = [latexstring(i-m-1) for i = 1:n]
	N = maximum(length.(labels))
	labels = lpad.(labels,N)
    for i=1:n
        f = Φ.ϕ[i]
        j = i-m-1
		ls = j > 0 ? ls=:dash : ls=:solid
        plot!(abs.(f.coefficients);label=labels[i],lw=3,ls=ls)
    end
    plot!(xlabel=L"n",ylabel=L"|c_n|",legend=:outertopright)
    plot!(;yscale,xlims,ylims,kwargs...)
end
