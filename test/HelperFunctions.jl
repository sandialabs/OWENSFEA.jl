using Test
import OWENSFEA

@testset "Gauss quadrature points" begin
    xi, weight = OWENSFEA.getGP(1)
    @test xi == [0.0]
    @test weight == [2.0]

    xi, weight = OWENSFEA.getGP(2)
    @test xi == [-sqrt(1 / 3), sqrt(1 / 3)]
    @test weight == [1.0, 1.0]

    xi, weight = OWENSFEA.getGP(3)
    @test xi == [-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)]
    @test weight == [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

    xi_inner = sqrt((3.0 - 2 * sqrt(6.0 / 5.0)) / 7.0)
    xi_outer = sqrt((3.0 + 2 * sqrt(6.0 / 5.0)) / 7.0)
    weight_inner = (18 + sqrt(30)) / 36.0
    weight_outer = (18 - sqrt(30)) / 36.0
    xi, weight = OWENSFEA.getGP(4)
    @test xi == [xi_inner, -xi_inner, xi_outer, -xi_outer]
    @test weight == [weight_inner, weight_inner, weight_outer, weight_outer]
end

@testset "Lagrange shape functions" begin
    N, p_N_x, Jac = OWENSFEA.calculateShapeFunctions(1, 0.25, [2.0, 6.0])
    @test N == [0.375, 0.625]
    @test p_N_x == [-0.25, 0.25]
    @test Jac == 2.0
    @test sum(N) == 1.0
    @test sum(p_N_x) == 0.0

    N, p_N_x, Jac = OWENSFEA.calculateShapeFunctions(2, 0.5, [2.0, 5.0, 8.0])
    @test N == [-0.125, 0.75, 0.375]
    @test isapprox(p_N_x, [0.0, -1.0 / 3.0, 1.0 / 3.0]; rtol=0.0, atol=eps())
    @test Jac == 3.0
    @test sum(N) == 1.0
    @test sum(p_N_x) == 0.0
end

@testset "Reduced DOF and boundary-condition maps" begin
    is_constrained = [0, 1, 0, 0, 1, 0]
    pBC = [1.0 2.0 10.0;
           3.0 1.0 -2.0]
    expected_map = [1.0, -1.0, 2.0, 3.0, -1.0, 4.0]

    @test OWENSFEA.calculateReducedDOFVector(3, 2, is_constrained) == [1, 3, 4, 6]
    @test OWENSFEA.calculateReducedDOFVector(1, 2, [1, 1]) == Int[]
    @test OWENSFEA.calculateBCMap(size(pBC, 1), pBC, 2, collect(1:6)) == expected_map
    @test OWENSFEA.constructReducedDispVectorMap(3, 2, 6, size(pBC, 1), pBC, is_constrained) == expected_map
end

@testset "Static boundary-condition application" begin
    K = [10.0 2.0 3.0 4.0;
         2.0 20.0 5.0 6.0;
         3.0 5.0 30.0 7.0;
         4.0 6.0 7.0 40.0]
    F = [1.0, 2.0, 3.0, 4.0]
    original_K = copy(K)
    original_F = copy(F)

    pBC = [1.0 2.0 1.5;
           2.0 1.0 -2.0]
    bc_map = [1.0, -1.0, -1.0, 2.0]
    BC = OWENSFEA.BC_struct(size(pBC, 1), pBC, 0, 0, [0.0, 1.0, 1.0, 0.0], bc_map, bc_map)

    K_bounded, F_bounded = OWENSFEA.applyBC(K, F, BC, 2)

    @test K_bounded == [10.0 0.0 0.0 4.0;
                        0.0 1.0 0.0 0.0;
                        0.0 0.0 1.0 0.0;
                        4.0 0.0 0.0 40.0]
    @test F_bounded == [4.0, 1.5, -2.0, 9.0]
    @test K == original_K
    @test F == original_F
end
