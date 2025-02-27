using Test

# Original:
# 100.116260 seconds (177.95 M allocations: 8.466 GiB, 3.35% gc time, 91.21% compilation time: <1% of which was recompilation)
# Tighter Mesh struct types
# 102.993400 seconds (177.25 M allocations: 8.409 GiB, 3.93% gc time, 96.15% compilation time: <1% of which was recompilation)
@time @testset "Cantilever Beam Modal" begin
    include("CantileverBeamModal.jl")
end

# # 50.008475 seconds (147.62 M allocations: 8.903 GiB, 5.53% gc time, 74.71% compilation time: <1% of which was recompilation)
# @time @testset "Cantilever Beam Deflection" begin
#     include("CantileverBeamDisplacement.jl")
# end

# # 73.630776 seconds (190.60 M allocations: 9.425 GiB, 6.37% gc time, 90.58% compilation time: <1% of which was recompilation)
# @time @testset "Cantilever Beam 45 degree sweep rotation" begin
#     include("CantileverBeamRotating.jl")
# end

# # 63.247943 seconds (169.29 M allocations: 9.411 GiB, 4.79% gc time, 86.57% compilation time: <1% of which was recompilation)
# @time @testset "Cantilever Beam MODAL 45 degree sweep with rotation" begin
#     include("CantileverBeamRotatingModal.jl")
# end

# # 297.712516 seconds (1.34 G allocations: 94.445 GiB, 5.81% gc time, 36.01% compilation time: <1% of which was recompilation)
# @time @testset "Cantilever Beam Unsteady Response to Unsteady Tip Load" begin
#     include("UnsteadyBeam.jl")
# end
