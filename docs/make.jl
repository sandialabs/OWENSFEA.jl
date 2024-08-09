using Documenter, Literate, OWENSFEA

# Build documentation
makedocs(;
    modules = [OWENSFEA],
    pages = [
        "Home" => "index.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSFEA.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)

deploydocs(
    repo = "github.com/sandialabs/OWENSFEA.jl.git",
)