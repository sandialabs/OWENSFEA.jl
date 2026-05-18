using Documenter, Literate, OWENSFEA

# Build documentation
makedocs(;
    modules = [OWENSFEA],
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Model Assembly" => "model_assembly.md",
        "Theory, Frames, and Units" => joinpath("theory", "frames_units.md"),
        "Validation and Testing" => "validation.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSFEA.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
)

deploydocs(
    repo = "github.com/sandialabs/OWENSFEA.jl.git",
)
