using TimerOutputs

reset_timer!()

@timeit "import packages" begin
    using GXBeam, LinearAlgebra
    import OWENSFEA
end


@timeit "setup" begin
    # include("../src/OWENSFEA.jl")
    include("./testdeps.jl")
    import Statistics
    using Test
    # Cantilever Beam/Rod (fixed - free)

    # Define Beam
    L = 0.5
    b = 0.05
    h = 0.02
    A = b * h
    E = 2.1e11
    v = 0.28
    G = E / 2(1 + v)
    rho = 7800.0
    Iyy = b * h^3 / 12
    Izz = b^3 * h / 12
    J = Iyy + Izz
    P = zeros(1)
    nev = 10 #modes
    ###############################################
    ######## Analytical
    ###############################################
    # Natural Frequencies
    c = sqrt(E * Iyy / (rho * A)) #for a beam

    Nmodes = 3
    k = [1.875, 4.694, 7.855]
    W = zeros(Nmodes)
    freqAnalytical = zeros(Nmodes)
    for n in 1:Nmodes
        W[n] = (2 * n + 1) * pi * c / (2 * L) # analytical natural frequency
        freqAnalytical[n] = sqrt(E * Iyy / (rho * A)) / (2 * pi) * (k[n] / L)^2
    end

    # println(W)
    # println(freqAnalytical)


    ###############################################
    ######## GXBeam
    ###############################################

    # create points
    # straight section of the beam
    r = [0.0, 0, 0]
    nelem = 40
    lengths, points, xm, Cab = discretize_beam(L, r, nelem)

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:(nelem + 1)

    # compliance matrix for each beam element
    compliance = fill(Diagonal([1 / (E * A), 1 / (G * A), 1 / (G * A), 1 / (G * J), 1 / (E * Iyy), 1 / (E * Izz)]), nelem)
    mass = fill(Diagonal([rho * A, rho * A, rho * A, rho * J, rho * Iyy, rho * Izz]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(
        points, start, stop;
        compliance = compliance,
        mass = mass,
        frames = Cab,
        lengths = lengths,
        midpoints = xm
    )

    system = GXBeam.DynamicSystem(assembly)
    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux = 0, uy = 0, uz = 0, theta_x = 0, theta_y = 0, theta_z = 0),
        # shear force on right tip
        nelem + 1 => PrescribedConditions(Fz_follower = 0.0)
    )

end

@timeit "initial compute" begin

    # eigenvalues and (right) eigenvectors
    # XXX This takes forever??
    @timeit "GXBeam.eigenvalue_analysis!" begin
        system, λ2, V, converged = eigenvalue_analysis!(
            system, assembly;
            prescribed_conditions = prescribed_conditions,
            nev = nev
        )
    end

    # corresponding left eigenvectors
    @timeit "GXBeam.left_eigenvectors" begin
        U = left_eigenvectors(system, λ2, V)
    end

    # post-multiply mass matrix with right eigenvector matrix
    # (we use this later for correlating eigenvalues)
    MV = system.M * V

    # process state and eigenstates
    @timeit "GXBeam.AssemblyState" begin
        state = AssemblyState(
            system, assembly;
            prescribed_conditions = prescribed_conditions
        )
        eigenstates = [
            AssemblyState(
                    V[:, k], system, assembly;
                    prescribed_conditions = prescribed_conditions
                ) for k in 1:nev
        ]
    end

    @timeit "GXBeam: postprocess" begin
        frequencyNative = [imag(λ2[k]) / (2 * pi) for k in 1:2:nev]

        # set previous left eigenvector matrix
        U_p = copy(U)

        # construct correlation matrix
        C = U_p * MV

        # correlate eigenmodes
        perm, corruption = correlate_eigenmodes(C)

        # re-arrange eigenvalues and eigenvectors
        λ2 = λ2[perm]
        U = U[perm, :]
        MV = MV[:, perm]
        eigenstates = eigenstates[perm]

        # update previous eigenvector matrix
        U_p .= U

        # update previous eigenvector matrix
        U_p .= U[1]

        freqGXBeam = [imag(λ2[k]) / (2 * pi) for k in 1:2:nev]
        # println(freqGXBeam)

        #These should match the FEA case from https://autofem.com/examples/determining_natural_frequencie.html
        # with frequencies of 67, 418, 1157
    end

end

#NOTE: with the full stiffness and mass matrices the GX modes matching the 2D analytical in plane modes are 1,3,5.

###############################################
######## OWENS
###############################################

@timeit "owens" begin
    # Create Mesh
    mesh, ort, joint = mesh_beam(;
        L1 = L / 2, #first section of beam length
        L2 = L / 2, #second section of beam length
        Nelem1 = 20,
        Nelem2 = 20,
        angleD = 0.0, # angle of second section of beam relative to first (0 is straight)
        zeroOffset = 0.0
    ) #r_b1[1]) #offset from 0 before first beam begins

    ort.Twist_d .= 180.0 #enforce flatwise orientation

    # Create Sectional Properties
    sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, length(mesh.z) - 1)

    for ii in 1:(length(mesh.z) - 1)
        ac = [0.0, 0.0]
        twist_d = [0.0, 0.0]
        rhoA = [rho * A, rho * A]
        EIyy = [E * Iyy, E * Iyy]
        EIzz = [E * Izz, E * Izz]
        GJ = [G * J, G * J]
        EA = [E * A, E * A]
        rhoIyy = [rho * Iyy, rho * Iyy]
        rhoIzz = [rho * Izz, rho * Izz]
        rhoJ = [rho * J, rho * J]
        zcm = [0.0, 0.0]
        ycm = [0.0, 0.0]
        a = [0.0, 0.0]

        EIyz = [0.0, 0.0]
        alpha1 = [0.0, 0.0] #This is always 0 in the element file, and it is unclear what it is used for since I can't find it being used in the code
        alpha2 = [0.0, 0.0]
        alpha3 = [0.0, 0.0]
        alpha4 = [0.0, 0.0]
        alpha5 = [0.0, 0.0]
        alpha6 = [0.0, 0.0]
        rhoIyz = [0.0, 0.0]
        local b = [0.0, 0.0]
        a0 = [0.0, 0.0]
        aeroCenterOffset = [0.0, 0.0]
        sectionPropsArray[ii] = OWENSFEA.SectionPropsArray(ac, twist_d, rhoA, EIyy, EIzz, GJ, EA, rhoIyy, rhoIzz, rhoJ, zcm, ycm, a, EIyz, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, rhoIyz, b, a0, aeroCenterOffset)
    end

    rotationalEffects = ones(mesh.numEl)
    el = OWENSFEA.El(sectionPropsArray, ort.Length, ort.Psi_d, ort.Theta_d, ort.Twist_d, rotationalEffects)

    # node, dof, boundary condition (bc)
    pBC = [
        1 1 0
        1 2 0
        1 3 0
        1 4 0
        1 5 0
        1 6 0
    ]

    feamodel = OWENSFEA.FEAModel(;
        analysisType = "M",
        dataOutputFilename = "none",
        joint = joint,
        gravityOn = false,
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        numNodes = mesh.numNodes
    )

    freq, damp, imagCompSign, U_x_0, U_y_0, U_z_0, theta_x_0, theta_y_0, theta_z_0, U_x_90, U_y_90, U_z_90, theta_x_90, theta_y_90, theta_z_90 = OWENSFEA.modal(feamodel, mesh, el; returnDynMatrices = true)

    # OWENS Frequencies that correspond to the GX beam are every other, and then 1,3,5 of the every other sets corresponds to the analytical

    # println("OWENS")
    # println(freq[1:10])
    freqOWENS = freq[1:2:10]
    freqOWENS2D = freqOWENS[1:2:end]
end

###############################################
######## TEST
###############################################

@timeit "TEST" begin
    for ifreq in 1:5
        atol = freqGXBeam[ifreq] * 0.01 # 1%
        @test isapprox(freqOWENS[ifreq], freqGXBeam[ifreq]; atol)
    end

    for ifreq in 1:3
        atol = freqAnalytical[ifreq] * 0.05 #5%
        @test isapprox(freqOWENS2D[ifreq], freqAnalytical[ifreq]; atol)
        # println(freqOWENS2D[ifreq])
        # println("error = $((freqOWENS2D[ifreq]-freqAnalytical[ifreq])/freqAnalytical[ifreq]*100)")
    end
end

print_timer()
