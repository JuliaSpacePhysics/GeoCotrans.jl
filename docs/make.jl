using Documenter
using GeoCotrans

DocMeta.setdocmeta!(GeoCotrans, :DocTestSetup, :(using GeoCotrans); recursive = true)

makedocs(
    sitename = "GeoCotrans.jl",
    format = Documenter.HTML(),
    modules = [GeoCotrans],
    pages = [
        "Home" => "index.md",
    ],
    checkdocs = :exports,
    doctest = true
)

deploydocs(
    repo = "github.com/juliaspacephysics/GeoCotrans.jl",
    push_preview = true
)
