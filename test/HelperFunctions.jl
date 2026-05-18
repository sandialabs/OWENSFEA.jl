using Test
using LinearAlgebra: I
import OWENSFEA

function thrown_message(f)
    try
        f()
    catch err
        return sprint(showerror, err)
    end
    return "NO_ERROR"
end

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

@testset "Element mass center offsets" begin
    rhoA = 10.0
    rhoIyy = 2.0
    rhoIzz = 3.0
    rhoIyz = 0.4
    rhoJ = 5.0
    ycm = 0.2
    zcm = -0.3
    x = 1.0
    y = 0.5
    z = 0.75
    integration_factor = 2.0

    M, Itens, xm = OWENSFEA.calculateElementMass(
        rhoA,
        rhoIyy,
        rhoIzz,
        rhoIyz,
        rhoJ,
        ycm,
        zcm,
        x,
        y,
        z,
        integration_factor,
        1.5,
        zeros(3, 3),
        zeros(3),
    )

    offset_y = y + ycm
    offset_z = z + zcm
    expected_mass = 21.5
    expected_inertia =
        rhoA * integration_factor *
        [
            offset_y^2+offset_z^2 -x*offset_y -x*offset_z
            -x*offset_y x^2+offset_z^2 -offset_y*offset_z
            -x*offset_z -offset_y*offset_z x^2+offset_y^2
        ] +
        integration_factor *
        [
            rhoJ 0.0 0.0
            0.0 rhoIyy rhoIyz
            0.0 rhoIyz rhoIzz
        ]

    @test M == expected_mass
    @test Itens == expected_inertia
    @test xm == rhoA * integration_factor * [x, offset_y, offset_z]
    @test xm[2] == 14.0
    @test xm[3] == 9.0
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

@testset "Joint dependent-active DOF maps" begin
    joint(joint_type; psi = 0.0, theta = 0.0, lx = 0.25, ly = 0.5, lz = 0.75) =
        [1.0, 1.0, 2.0, joint_type, lx, ly, lz, theta]

    expected_by_type = Dict(
        0 => (d = [7, 8, 9, 10, 11, 12], a = [1, 2, 3, 4, 5, 6],
              T = Matrix{Float64}(I, 6, 6)),
        1 => (d = [7, 8, 9], a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 11.0, 12.0],
              T = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0]),
        2 => (d = [7, 8, 9, 10, 12], a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 11.0],
              T = [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 1.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 0.0 1.0 0.0]),
        3 => (d = [7, 8, 9, 11, 12], a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0],
              T = [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 0.0 1.0 0.0]),
        4 => (d = [7, 8, 9, 10, 11], a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 12.0],
              T = [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 1.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 1.0 0.0 0.0]),
        5 => (d = [7, 8, 9, 10, 11, 12], a = [1, 2, 3, 4, 5, 6],
              T = [1.0 0.0 0.0 0.0 0.75 -0.5;
                   0.0 1.0 0.0 -0.75 0.0 0.25;
                   0.0 0.0 1.0 0.5 -0.25 0.0;
                   0.0 0.0 0.0 1.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 1.0 0.0;
                   0.0 0.0 0.0 0.0 0.0 1.0]),
    )

    for joint_type in 0:5
        Tda, dDOF, aDOF = OWENSFEA.createTda(joint_type, 2, 1, 0.0, 0.0, joint(joint_type))
        expected = expected_by_type[joint_type]
        @test dDOF == expected.d
        @test vec(aDOF) == expected.a
        @test Matrix(Tda) == expected.T
        @test eltype(Matrix(Tda)) === Float64
    end

    Tda, dDOF, aDOF = OWENSFEA.createTda(4, 2, 1, 0.0, 90.0, joint(4; theta = 90.0))
    @test dDOF == [7, 8, 9, 11, 12]
    @test vec(aDOF) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0]
    @test Matrix(Tda) == expected_by_type[3].T

    @test thrown_message(() -> OWENSFEA.createTda(9, 2, 1, 0.0, 0.0, joint(9))) ==
        "Correct jointType not specified, should be 1, 2, 3, 4, or 5"
    @test thrown_message(() -> OWENSFEA.getNodeMaps(zeros(2, 2), zeros(2, 3), 1, 2, [1, 2], [1, 2, 3], Int[])) ==
        "Singular joint transformation matrix. Exiting"
end

@testset "Joint slave DOF extraction branches" begin
    jointrow(joint_type; psi = 0.0, theta = 0.0) = [1.0 1.0 2.0 joint_type 0.0 0.0 psi theta]
    expected_slave_dofs = (
        (jointrow(1), [7, 8, 9]),
        (jointrow(2; psi = 90.0), [7, 8, 9, 11, 12]),
        (jointrow(2; psi = 0.0), [7, 8, 9, 10, 12]),
        (jointrow(3; theta = 90.0), [7, 8, 9, 10, 11]),
        (jointrow(3; psi = 90.0), [7, 8, 9, 10, 12]),
        (jointrow(3), [7, 8, 9, 11, 12]),
        (jointrow(4; theta = 90.0, psi = 90.0), [7, 8, 9, 10, 12]),
        (jointrow(4; theta = 90.0), [7, 8, 9, 11, 12]),
        (jointrow(4), [7, 8, 9, 10, 11]),
        (jointrow(5), [7, 8, 9, 10, 11, 12]),
    )

    for (joint_data, slave_dofs) in expected_slave_dofs
        adNumDof, aNumDof, actual_slave_dofs = OWENSFEA.extractdaInfo(joint_data, 2, 6)
        @test adNumDof == 12
        @test aNumDof == 12 - length(slave_dofs)
        @test actual_slave_dofs == slave_dofs
    end
end

@testset "Concentrated nodal term parsing" begin
    full_terms = Any[
        1 "M6" 2 3 4.5;
        1 "K6" 3 4 5.5;
        2 "C6" 4 5 6.5;
        2 "F" 6 1 7.5
    ]
    nodal = @test_logs (:warn, "Concentrated loads are one-dimensional, second degree of freedom input is ignored.") OWENSFEA.applyConcentratedTerms(2, 6; data = full_terms)
    @test nodal.concMass[2, 3] == 4.5
    @test nodal.concStiff[3, 4] == 5.5
    @test nodal.concDamp[10, 11] == 6.5
    @test nodal.concLoad[12] == 7.5
    @test size(nodal.concMass) == (12, 12)
    @test eltype(nodal.concLoad) === Float64

    diagonal_terms = Any[
        1 "M" 2 1.5;
        1 "K" 3 2.5;
        1 "C" 4 3.5;
        1 "F" 5 4.5
    ]
    diagonal = @test_logs (:warn, "Only one degree of freedom given for concentrated mass; applying at diagonal.") (:warn, "Only one degree of freedom given for concentrated stiffness; applying at diagonal.") (:warn, "Only one degree of freedom given for concentrated damping; applying at diagonal.") OWENSFEA.applyConcentratedTerms(1, 6; data = diagonal_terms)
    @test diagonal.concMass[2, 2] == 1.5
    @test diagonal.concStiff[3, 3] == 2.5
    @test diagonal.concDamp[4, 4] == 3.5
    @test diagonal.concLoad[5] == 4.5
    @test thrown_message(() -> OWENSFEA.applyConcentratedTerms(1, 6; data = Any[1 "X" 1 1.0])) ==
        "Unknown Nodal Data Type"
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
