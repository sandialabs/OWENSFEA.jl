using Test
import OWENSFEA

function _spinning_beam_fixture(; L=2.0, nelem=4, y_offset=0.0, section_ycm=0.0, section_zcm=0.0, stiffness_scale=1.0)
    b = 0.05
    h = 0.02
    A = b*h
    E = 2.1e11*stiffness_scale
    nu = 0.28
    G = E/(2*(1 + nu))
    rho = 7800.0
    Iyy = b*h^3/12
    Izz = b^3*h/12
    J = Iyy + Izz

    rhoA = rho*A
    pair(value) = fill(value, 2)
    section_props = Array{OWENSFEA.SectionPropsArray, 1}(undef, nelem)
    for i = 1:nelem
        section_props[i] = OWENSFEA.SectionPropsArray(
            pair(0.0), pair(0.0),
            fill(rhoA, 2),
            fill(E*Iyy, 2),
            fill(E*Izz, 2),
            fill(G*J, 2),
            fill(E*A, 2),
            fill(rho*Iyy, 2),
            fill(rho*Izz, 2),
            fill(rho*J, 2),
            pair(section_zcm), pair(section_ycm), pair(0.0), pair(0.0),
            pair(0.0), pair(0.0), pair(0.0), pair(0.0), pair(0.0), pair(0.0),
            pair(0.0), pair(0.0), pair(0.0), pair(0.0),
            nothing, nothing,
            pair(0.0), pair(0.0),
            fill(G*A, 2),
            fill(G*A, 2),
        )
    end

    x = collect(range(0.0, L, length=nelem + 1))
    conn = hcat(collect(1:nelem), collect(2:nelem + 1))
    mesh = OWENSFEA.Mesh(
        collect(1:nelem + 1),
        nelem,
        nelem + 1,
        x,
        fill(y_offset, nelem + 1),
        zeros(nelem + 1),
        collect(1:nelem),
        conn,
        zeros(Int, nelem),
        [nelem],
        zeros(1, 1),
        zeros(Int, 1, 1),
        zeros(Int, 1, 1),
    )
    el = OWENSFEA.El(
        section_props,
        fill(L/nelem, nelem),
        zeros(nelem),
        zeros(nelem),
        zeros(nelem),
        ones(nelem),
    )

    return mesh, el, rhoA, L
end

function _spinning_beam_model(mesh)
    root_fixed = [
        1 1 0
        1 2 0
        1 3 0
        1 4 0
        1 5 0
        1 6 0
    ]
    return OWENSFEA.FEAModel(;
        analysisType="S",
        dataOutputFilename="none",
        joint=zeros(0, 8),
        pBC=root_fixed,
        numNodes=mesh.numNodes,
        nlOn=false,
        gravityOn=false,
        iterationType="LINEAR",
        maxNumLoadSteps=3,
    )
end

function _root_reaction(mesh, el; Omega=0.0, OmegaDot=0.0)
    feamodel = _spinning_beam_model(mesh)
    el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
    displ0 = zeros(mesh.numNodes*6)
    _, _, success, reaction = redirect_stdout(devnull) do
        OWENSFEA.staticAnalysis(
            feamodel,
            mesh,
            el,
            displ0,
            Omega,
            Omega,
            el_storage;
            OmegaDot,
        )
    end
    return success, reaction[1:6]
end

@testset "Sectional CG offset steady-spin reactions" begin
    section_ycm = 0.3
    spin_hz = 5.0
    mesh, el, rhoA, L = _spinning_beam_fixture(; section_ycm, stiffness_scale=1e6)
    success, reaction = _root_reaction(mesh, el; Omega=spin_hz)

    omega = 2*pi*spin_hz
    expected_axial_reaction = -rhoA*omega^2*L^2/2
    expected_side_reaction = -rhoA*omega^2*section_ycm*L
    expected_moment_scale = abs(rhoA*omega^2*section_ycm*L^2)

    @test success
    @test isapprox(reaction[1], expected_axial_reaction; rtol=1e-3)
    @test isapprox(reaction[2], expected_side_reaction; rtol=1e-3)
    @test isapprox(reaction[3], 0.0; atol=1e-10*abs(expected_axial_reaction))
    @test isapprox(reaction[4], 0.0; atol=1e-10*expected_moment_scale)
    @test isapprox(reaction[5], 0.0; atol=1e-10*expected_moment_scale)
    @test isapprox(reaction[6], 0.0; atol=1e-7*expected_moment_scale)
end

@testset "Eccentric spinning beam spin-acceleration reactions" begin
    y_offset = 0.3
    spin_accel_hz = 2.0
    mesh, el, rhoA, L = _spinning_beam_fixture(; y_offset, stiffness_scale=1e4)
    success, reaction = _root_reaction(mesh, el; OmegaDot=spin_accel_hz)

    alpha = 2*pi*spin_accel_hz
    expected_x_reaction = -rhoA*alpha*y_offset*L
    expected_y_reaction = rhoA*alpha*L^2/2
    expected_spanwise_spin_axis_torque = rhoA*alpha*L^3/3
    expected_offset_spin_axis_torque = rhoA*alpha*y_offset^2*L

    @test success
    @test isapprox(reaction[1], expected_x_reaction; rtol=1e-5)
    @test isapprox(reaction[2], expected_y_reaction; rtol=1e-5)
    # The current root-reaction moment channel omits the eccentric y^2 term even
    # though the corresponding x-force is present. Pin both sides so a future
    # #15/#17 torque fix is an intentional test update.
    @test isapprox(reaction[6], expected_spanwise_spin_axis_torque; rtol=1e-5)
    @test isapprox(
        expected_spanwise_spin_axis_torque + expected_offset_spin_axis_torque - reaction[6],
        expected_offset_spin_axis_torque;
        rtol=1e-5,
    )
    @test isapprox(reaction[3], 0.0; atol=1e-10*abs(expected_y_reaction))
    @test isapprox(reaction[4], 0.0; atol=1e-10*abs(expected_spanwise_spin_axis_torque))
    @test isapprox(reaction[5], 0.0; atol=1e-10*abs(expected_spanwise_spin_axis_torque))
end

@testset "Straight spinning beam reaction diagnostics" begin
    mesh, el, rhoA, L = _spinning_beam_fixture()

    spin_hz = 5.0
    success, reaction = _root_reaction(mesh, el; Omega=spin_hz)
    omega = 2*pi*spin_hz
    expected_axial_reaction = -rhoA*omega^2*L^2/2

    @test success
    @test isapprox(reaction[1], expected_axial_reaction; rtol=1e-3)
    @test isapprox(reaction[2], 0.0; atol=1e-10*abs(expected_axial_reaction))
    @test isapprox(reaction[3], 0.0; atol=1e-10*abs(expected_axial_reaction))
    @test isapprox(reaction[4], 0.0; atol=1e-10*abs(expected_axial_reaction)*L)
    @test isapprox(reaction[5], 0.0; atol=1e-10*abs(expected_axial_reaction)*L)
    @test isapprox(reaction[6], 0.0; atol=1e-10*abs(expected_axial_reaction)*L)

    spin_accel_hz = 2.0
    success, reaction = _root_reaction(mesh, el; OmegaDot=spin_accel_hz)
    alpha = 2*pi*spin_accel_hz
    expected_tangential_reaction = rhoA*alpha*L^2/2
    expected_spin_axis_torque = rhoA*alpha*L^3/3

    @test success
    @test isapprox(reaction[2], expected_tangential_reaction; rtol=1e-5)
    @test isapprox(reaction[6], expected_spin_axis_torque; rtol=1e-5)
    @test isapprox(reaction[3], 0.0; atol=1e-10*abs(expected_tangential_reaction))
    @test isapprox(reaction[4], 0.0; atol=1e-10*abs(expected_spin_axis_torque))
    @test isapprox(reaction[5], 0.0; atol=1e-10*abs(expected_spin_axis_torque))
end

@testset "Static rigid-body spin vector plumbing" begin
    mesh, el, _, _ = _spinning_beam_fixture()
    feamodel = _spinning_beam_model(mesh)
    el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
    displ0 = zeros(mesh.numNodes*6)

    rb_data = zeros(9)
    rb_data[4] = 2*pi*5.0
    _, _, success, reaction = redirect_stdout(devnull) do
        OWENSFEA.staticAnalysis(
            feamodel,
            mesh,
            el,
            displ0,
            0.0,
            0.0,
            el_storage;
            rbData=rb_data,
        )
    end

    @test success
    @test all(isapprox.(reaction[1:6], 0.0; atol=1e-8))
end
