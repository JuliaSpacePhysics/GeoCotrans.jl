using Documenter
using GeoCotrans

DocMeta.setdocmeta!(GeoCotrans, :DocTestSetup, :(using GeoCotrans); recursive = true)

makedocs(
    sitename = "GeoCotrans.jl",
    format = Documenter.HTML(),
    modules = [GeoCotrans],
    pages = [
        "Home" => "index.md",
        "Validation & Comparison" => "coords.md"
    ],
    checkdocs = :exports,
    doctest = true
)

deploydocs(
    repo = "github.com/JuliaSpacePhysics/GeoCotrans.jl",
    push_preview = true
)
