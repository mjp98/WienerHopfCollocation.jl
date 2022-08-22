using WienerHopfCollocation
using Documenter

DocMeta.setdocmeta!(WienerHopfCollocation, :DocTestSetup, :(using WienerHopfCollocation); recursive=true)

makedocs(;
    modules=[WienerHopfCollocation],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/WienerHopfCollocation.jl/blob/{commit}{path}#{line}",
    sitename="WienerHopfCollocation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/WienerHopfCollocation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/WienerHopfCollocation.jl",
    devbranch="main",
)
