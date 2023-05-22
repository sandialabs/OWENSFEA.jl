using Documenter, Literate, GyricFEA

# Build documentation
makedocs(;
    modules = [GyricFEA],
    pages = [
        "Home" => "index.md",
        # "Getting Started" => joinpath("examples", "guide.md"),
        # "Section Properties and Strain Recovery" => joinpath("examples", "section.md"),
        # "Sensitivity Analysis" => joinpath("examples", "sensitivities.md"),
        # "DifferentialEquations" => joinpath("examples", "diffeq.md"),
        # "Examples" => [
        #     joinpath("examples", "cantilever.md"),
        #     joinpath("examples", "overdetermined.md"),
        #     joinpath("examples", "tipforce.md"),
        #     joinpath("examples", "tipmoment.md"),
        #     joinpath("examples", "curved.md"),
        #     joinpath("examples", "rotating.md"),
        #     joinpath("examples", "excited.md"),
        #     joinpath("examples", "wind-turbine-blade.md"),
        #     joinpath("examples", "static-joined-wing.md"),
        #     joinpath("examples", "dynamic-joined-wing.md"),
        #     joinpath("examples", "vertical-axis-wind-turbine.md"),
        # ],
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "GyricFEA.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)