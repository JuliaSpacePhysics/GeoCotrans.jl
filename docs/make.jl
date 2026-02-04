using Documenter
using GeoCotrans

DocMeta.setdocmeta!(GeoCotrans, :DocTestSetup, :(using GeoCotrans; using GeoCotrans.FieldLineTracing); recursive = true)

makedocs(
    sitename = "GeoCotrans.jl",
    format = Documenter.HTML(),
    modules = [GeoCotrans, GeoCotrans.FieldLineTracing],
    pages = [
        "Home" => "index.md",
        "Validation & Comparison" => "coords.md",
        "Field Line Tracing" => "field_line_tracing.md"
    ],
    checkdocs = :exports,
    doctest = true
)

deploydocs(
    repo = "github.com/JuliaSpacePhysics/GeoCotrans.jl",
    push_preview = true
)
