using Test

@testset "Cantilever Beam Modal" begin
    include("CantileverBeamModal.jl")
end

@testset "Cantilever Beam Deflection" begin
    include("CantileverBeamDisplacement.jl")
end

@testset "Cantilever Beam 45 degree sweep rotation" begin
    include("CantileverBeamRotating.jl")
end

@testset "Cantilever Beam MODAL 45 degree sweep with rotation" begin
    include("CantileverBeamRotatingModal.jl")
end

@testset "Cantilever Beam Unsteady Response to Unsteady Tip Load" begin
    include("UnsteadyBeam.jl")
end
