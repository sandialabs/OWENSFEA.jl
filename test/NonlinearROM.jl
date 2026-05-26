using LinearAlgebra
using Test

import OWENSFEA

function straight_cantilever_fixture(; nelem=8, length=2.0)
    x = collect(range(0.0, length; length=nelem + 1))
    y = zeros(nelem + 1)
    z = zeros(nelem + 1)
    conn = hcat(collect(1:nelem), collect(2:nelem + 1))

    mesh = OWENSFEA.Mesh(
        collect(1:nelem + 1),
        nelem,
        nelem + 1,
        x,
        y,
        z,
        collect(1:nelem),
        conn,
        zeros(Int, nelem),
        [nelem],
        zeros(1, 1),
        zeros(Int, 1, 1),
        zeros(Int, 1, 1),
    )

    width = 0.045
    height = 0.018
    area = width * height
    youngs_modulus = 70.0e9
    poisson = 0.33
    shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson))
    density = 2700.0
    iyy = width * height^3 / 12.0
    izz = width^3 * height / 12.0
    polar_i = iyy + izz

    props = Vector{OWENSFEA.SectionPropsArray}(undef, nelem)
    for i in eachindex(props)
        ac = [0.0, 0.0]
        twist = [0.0, 0.0]
        rhoA = [density * area, density * area]
        EIyy = [youngs_modulus * iyy, youngs_modulus * iyy]
        EIzz = [youngs_modulus * izz, youngs_modulus * izz]
        GJ = [shear_modulus * polar_i, shear_modulus * polar_i]
        EA = [youngs_modulus * area, youngs_modulus * area]
        rhoIyy = [density * iyy, density * iyy]
        rhoIzz = [density * izz, density * izz]
        rhoJ = [density * polar_i, density * polar_i]
        zero2 = [0.0, 0.0]
        GAy = [5.0 / 6.0 * shear_modulus * area, 5.0 / 6.0 * shear_modulus * area]
        GAz = copy(GAy)
        props[i] = OWENSFEA.SectionPropsArray(
            ac,
            twist,
            rhoA,
            EIyy,
            EIzz,
            GJ,
            EA,
            rhoIyy,
            rhoIzz,
            rhoJ,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            zero2,
            nothing,
            nothing,
            zero2,
            zero2,
            GAy,
            GAz,
        )
    end

    el = OWENSFEA.El(
        props,
        fill(length / nelem, nelem),
        zeros(nelem),
        zeros(nelem),
        zeros(nelem),
        ones(nelem),
    )

    return mesh, el, zeros(0, 8)
end

function cantilever_feamodel(analysis_type, mesh, joint; nl_on=true, num_modes=36)
    pBC = [
        1 1 0
        1 2 0
        1 3 0
        1 4 0
        1 5 0
        1 6 0
    ]

    return OWENSFEA.FEAModel(;
        analysisType=analysis_type,
        dataOutputFilename="none",
        joint,
        pBC,
        numNodes=mesh.numNodes,
        nlOn=nl_on,
        gravityOn=false,
        iterationType="DI",
        tolerance=1.0e-9,
        maxIterations=100,
        RayleighAlpha=0.0,
        RayleighBeta=0.0,
        numModes=num_modes,
    )
end

function run_transient_fixture(analysis_type; nl_on=true)
    mesh, el, joint = straight_cantilever_fixture()
    feamodel = cantilever_feamodel(analysis_type, mesh, joint; nl_on)
    el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
    ndof = mesh.numNodes * 6
    rom = analysis_type == "ROM" ? OWENSFEA.reducedOrderModel(el_storage, feamodel, mesh, el, zeros(ndof)) : nothing
    disp_data = OWENSFEA.DispData(zeros(ndof), zeros(ndof), zeros(ndof), zeros(ndof))

    delta_t = 0.002
    times = delta_t:delta_t:0.04
    tip_z_dof = (mesh.numNodes - 1) * 6 + 3
    tip_history = zeros(length(times))

    for (itime, time) in enumerate(times)
        load = 200.0 * sin(pi * time / last(times))
        external_force = [load]
        force_dof = [tip_z_dof]
        if analysis_type == "ROM"
            _, disp_out, _ = OWENSFEA.structuralDynamicsTransientROM(
                feamodel,
                mesh,
                el,
                disp_data,
                0.0,
                0.0,
                time,
                delta_t,
                el_storage,
                rom,
                external_force,
                force_dof,
                I(3),
                zeros(9),
            )
        else
            _, disp_out, _ = OWENSFEA.structuralDynamicsTransient(
                feamodel,
                mesh,
                el,
                disp_data,
                0.0,
                0.0,
                time,
                delta_t,
                el_storage,
                external_force,
                force_dof,
                I(3),
                zeros(9),
            )
        end

        if analysis_type == "ROM"
            disp_data = OWENSFEA.DispData(
                disp_out.displ_sp1,
                disp_out.displdot_sp1,
                disp_out.displddot_sp1,
                disp_data.displ_s,
                disp_out.eta_sp1,
                disp_out.etadot_sp1,
                disp_out.etaddot_sp1,
            )
        else
            disp_data = OWENSFEA.DispData(
                disp_out.displ_sp1,
                disp_out.displdot_sp1,
                disp_out.displddot_sp1,
                disp_data.displ_s,
            )
        end
        tip_history[itime] = disp_out.displ_sp1[tip_z_dof]
    end

    return tip_history
end

full_linear = run_transient_fixture("TNB"; nl_on=false)
full_nonlinear = run_transient_fixture("TNB"; nl_on=true)
rom_nonlinear = run_transient_fixture("ROM"; nl_on=true)

@test maximum(abs.(full_nonlinear)) > 1.0e-3
@test norm(full_nonlinear - full_linear) / norm(full_linear) > 1.0e-4
@test norm(rom_nonlinear - full_nonlinear) / norm(full_nonlinear) < 5.0e-2
@test isapprox(rom_nonlinear[end], full_nonlinear[end]; rtol=5.0e-2, atol=5.0e-5)
