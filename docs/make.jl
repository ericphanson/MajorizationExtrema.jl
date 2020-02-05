using Documenter, MajorizationExtrema

makedocs(;
    modules=[MajorizationExtrema],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        # "Subentropy counterexample" => "subentropy.md",
    ],
    repo="https://github.com/ericphanson/MajorizationExtrema.jl/blob/{commit}{path}#L{line}",
    sitename="MajorizationExtrema.jl",
    authors="Eric",
    assets=String[],
)

deploydocs(;
    repo="github.com/ericphanson/MajorizationExtrema.jl",
)
