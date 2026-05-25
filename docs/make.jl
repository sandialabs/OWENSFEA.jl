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
    format = Documenter.HTML(
        repolink = "https://github.com/sandialabs/OWENSFEA.jl",
        edit_link = "master",
    ),
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/sandialabs/OWENSFEA.jl.git",
        devbranch = "master",
    )
end
