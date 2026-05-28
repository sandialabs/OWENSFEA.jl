using Documenter, Literate, OWENSFEA

# Build documentation
makedocs(;
    modules = [OWENSFEA],
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Test-backed Examples" => "examples.md",
        "Model Assembly" => "model_assembly.md",
        "Theory, Frames, and Units" => joinpath("theory", "frames_units.md"),
        "Validation and Testing" => "validation.md",
        "Developer Guide" => "developer_guide.md",
        "Reference" => [
            "API Map" => joinpath("reference", "reference.md"),
            "Autodocs by Source" => joinpath("reference", "autodocs.md"),
        ],
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
