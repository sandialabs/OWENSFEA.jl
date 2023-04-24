using GXBeam, LinearAlgebra
using Test
include("./testdeps.jl")
import Statistics
import GyricFEA
# include("../src/GyricFEA.jl")
# Cantilever Beam/Rod (fixed - free)

function runme(islinear)

    # Define Beam
    L = 0.5
    b = 0.05
    h = 0.02
    A = b*h
    E = 2.1e11
    v = 0.28
    G = E/2(1+v)
    rho = 7800.0
    Iyy = b*h^3/12
    Izz = b^3*h/12
    J = Iyy+Izz
    P = zeros(1)
    nev = 10 #modes
    sweepD = 45.0
    ###############################################
    ######## Analytical
    ###############################################
    RPM = 0:1000:6000
    dispAnalytical = zeros(length(RPM))

    for i = 1:length(RPM)
        dispAnalytical[i] = RPM[i]*L^3/(3*E*Iyy)
    end
    # println(dispAnalytical)
    ###############################################
    ######## GXBeam
    ###############################################

    # create points
    # straight section of the beam
    # straight section of the beam
    L_b1 = L/2 # inch
    r_b1 = [0.0, 0, 0]
    nelem_b1 = 20
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # swept section of the beam
    L_b2 = L/2 # inch
    r_b2 = [L/2, 0, 0]
    nelem_b2 = 20
    cs, ss = cos(sweepD*pi/180), sin(sweepD*pi/180)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2;
        frame = frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem_b1 + nelem_b2
    stop = 2:nelem_b1 + nelem_b2 + 1
    lengths = vcat(lengths_b1, lengths_b2)
    xm = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    # compliance matrix for each beam element
    compliance = fill(Diagonal([1/(E*A), 1/(G*A), 1/(G*A), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)
    mass = fill(Diagonal([rho*A, rho*A, rho*A, rho*J, rho*Iyy, rho*Izz]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop;
    compliance = compliance,
    mass = mass,
    frames = Cab,
    lengths = lengths,
    midpoints = xm)

    # run an analysis for each prescribed tip load
    states = Vector{AssemblyState{Float64}}(undef, length(RPM))
    λ2 = Vector{Vector{ComplexF64}}(undef,length(RPM))
    U = Vector{Matrix{ComplexF64}}(undef, length(RPM))
    MV = Vector{Matrix{ComplexF64}}(undef, length(RPM))
    state = Vector{AssemblyState{Float64}}(undef, length(RPM))
    eigenstates = Vector{Vector{AssemblyState{ComplexF64}}}(undef, length(RPM))
    deformedxyz = zeros(length(points),3,length(RPM))

    for i = 1:length(RPM)
        system = GXBeam.StaticSystem(assembly)
        w0 = [0, 0, RPM[i]*(2*pi)/60]
        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

        # perform a static analysis
        system, converged = steady_state_analysis(assembly,
            angular_velocity = w0,
            prescribed_conditions = prescribed_conditions,
            linear = islinear)

        # process state and eigenstates
        states[i] = AssemblyState(system, assembly;
        prescribed_conditions = prescribed_conditions)

        for (ii,point) in enumerate(states[i].points)
            deformedxyz[ii,:,i] = point.u
        end
    end

    # for ilevel = 1:length(RPM)
    #     println(deformedxyz[end,3,ilevel])
    # end

    ###############################################
    ######## OWENS
    ###############################################

    # Create Mesh
    mesh, ort, joint = mesh_beam(;L1 = L/2, #first section of beam length
    L2 = L/2, #second section of beam length
    Nelem1 = 20,
    Nelem2 = 20,
    angleD = sweepD, # angle of second section of beam relative to first (0 is straight)
    zeroOffset = 0.0,
    vertical = false)#r_b1[1]) #offset from 0 before first beam begins

    # Create Sectional Properties
    sectionPropsArray = Array{GyricFEA.SectionPropsArray, 1}(undef, length(mesh.z)-1)

    for ii = 1:length(mesh.z)-1
        ac = [0.0,0.0]
        twist_d = [0.0,0.0]
        rhoA = [rho*A,rho*A]
        EIzz = [E*Iyy,E*Iyy]
        EIyy = [E*Izz,E*Izz]
        GJ = [G*J,G*J]
        EA = [E*A,E*A]
        rhoIyy = [rho*Iyy,rho*Iyy]
        rhoIzz = [rho*Izz,rho*Izz]
        rhoJ = [rho*J,rho*J]
        zcm = [0.0,0.0]
        ycm = [0.0,0.0]
        a = [0.0,0.0]

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
        sectionPropsArray[ii] = GyricFEA.SectionPropsArray(ac,twist_d,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)
    end

    rotationalEffects = ones(mesh.numEl)
    el = GyricFEA.El(sectionPropsArray,ort.Length,ort.Psi_d,ort.Theta_d,ort.Twist_d,rotationalEffects)

    M25_beam = zeros(length(RPM))
    Mcurv_beam = zeros(length(RPM))
    Msweep_beam = zeros(length(RPM))
    Fn_beam = zeros(length(RPM))
    Ft_beam = zeros(length(RPM))
    Fz_beam = zeros(length(RPM))
    Ux_beam = zeros(length(RPM))
    Uy_beam = zeros(length(RPM))
    Uz_beam = zeros(length(RPM))
    Θ25_beam = zeros(length(RPM))
    Θcurv_beam = zeros(length(RPM))
    Θsweep_beam = zeros(length(RPM))


    # node, dof, boundary condition (bc)
    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0]

    numDOFPerNode = 6
    delta_t = 0.1
    finalt = 0.1
    numTS = round(Int,finalt/delta_t)
    TOL = 1e-3
    uHist = zeros(length(RPM),mesh.numNodes*numDOFPerNode,numTS+1) #put here so it keeps the last RRPMM solution in the scope


    for (ii,rpm) in enumerate(RPM)

        feamodel = GyricFEA.FEAModel(;analysisType = "TNB",
        outFilename = "none",
        joint,
        pBC = pBC,
        numNodes = mesh.numNodes,
        nlOn = !islinear,
        tolerance = 1.0e-06,
        maxIterations = 500,
        maxNumLoadSteps = 20,
        gravityOn = false,
        iterationType="DI",
        aeroElasticOn=false)

        elStorage = GyricFEA.initialElementCalculations(feamodel,el,mesh)

        structureMass,structureMOI,structureMassCenter = GyricFEA.calculateStructureMassProps(elStorage)

        u_s = zeros(mesh.numNodes*numDOFPerNode,1)
        # u_s[:] = 2.0.*[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0002786498297864991, 0.0001929465203607998, 0.0, 0.0, 0.0, 0.0, 0.0005529717288725646, 0.0009313743703659432, 0.0, 0.0, 0.0, 0.0, 0.0008206037989981936, 0.0024888842441722087, 0.0, 0.0, 0.0, 0.0, 0.001078899807277038, 0.005145350307976802, 0.0, 0.0, 0.0, 0.0, 0.001324819334136119, 0.009188670895859945, 0.0, 0.0, 0.0, 0.0, 0.0015548051834367155, 0.014916237961941218, 0.0, 0.0, 0.0, 0.0, 0.0017646462583511772, 0.022636039240240186, 0.0, 0.0, 0.0, 0.0, 0.0019493252555308333, 0.03266729739051467, 0.0, 0.0, 0.0, 0.0, 0.0021028521403319293, 0.04534054000670048, 0.0, 0.0, 0.0, 0.0, 0.002218086539899789, 0.060996983993136356, 0.0, 0.0, 0.0, 0.0, 0.0022865550199244592, 0.07998710843940031, 0.0, 0.0, 0.0, 0.0, 0.00229827275383148, 0.10266828294762664, 0.0, 0.0, 0.0, 0.0, 0.0022415833476843727, 0.12940131483878864, 0.0, 0.0, 0.0, 0.0, 0.0022415833476843727, 0.12940131483878864, 0.0, 0.0, 0.0, 0.0, 0.019852328157362865, 0.14717629328985662, 0.0, 0.0, 0.0, 0.0, 0.038470151060156624, 0.16600657787773634, 0.0, 0.0, 0.0, 0.0, 0.057394650935991415, 0.1851755162598843, 0.0, 0.0, 0.0, 0.0]
        # u_s = GyricFEA.setInitialConditions(feamodel.initCond,u_s,numDOFPerNode)
        udot_s = u_s*0
        uddot_s = u_s*0

        u_j=u_s
        udot_j = 0.0
        uddot_j = 0.0
        azi_s = 0
        azi_j = azi_s

        Omega_s = rpm/60.0#*2*pi
        Omega_j = Omega_s
        OmegaDot_j = 0.0
        Fexternal = []
        Fdof = []

        elStrain = fill(GyricFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)
        dispOut = GyricFEA.DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492)) #TODO: not hard coded
        FReactionHist = zeros(numTS+1,6)
        aziHist = zeros(numTS+1)

        FReactionsm1 = zeros(6)
        FReactionHist[1,:] = FReactionsm1
        FReaction_j = FReactionsm1

        # Run unsteady analysis

        for (i,time) in enumerate(delta_t:delta_t:finalt)
            # println("Time: $time")

            iwhile = 0
            uNorm = 1e5
            aziNorm = 1e5
            while uNorm > TOL || aziNorm > TOL #TODO: add TOL to feamodel
                # println("While Loop Iteration: $iwhile")

                iwhile += 1
                dispData = GyricFEA.DispData(u_s,udot_s,uddot_s,0.0)


                CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]
                CN2P=1.0*LinearAlgebra.I(3)
                CN2H = CP2H*CN2P

                rbData = zeros(9) #TODO

                #TODO: time isn't used in the function???
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,time,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)

                ## update timestepping variables and other states, store in history arrays
                u_jLast = u_j
                u_j = dispOut.displ_sp1
                udot_j  = dispOut.displdot_sp1
                uddot_j = dispOut.displddot_sp1

                azi_jLast = azi_j
                azi_j = azi_s + Omega_j*delta_t*2*pi

                uNorm = LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)
                aziNorm = LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)
                println("uNorm: $uNorm axiNorm: $aziNorm")
            end

            u_s = u_j
            udot_s = udot_j
            uddot_s = uddot_j
            azi_s = azi_j
            Omega_s = Omega_j

            uHist[ii,:,i+1] = u_s
            FReactionHist[i+1,:] = FReaction_j
            aziHist[i+1] = azi_s

        end # time stepping

        # Rotating Frame of Reference, 6 DOF where
        # 1 = turbine vertical force
        # 2 = turbine 2D slice tangential force
        # 3 = turbine 2D slice normal force
        # 4 = blade M25 twisting moment
        # 5 = blade curvature twisting moment  or is this sweep?
        # 6 = blade sweep moment.  and this is curvature moment?

        Fn_beam[ii] = FReactionHist[end,3]
        Ft_beam[ii] = FReactionHist[end,2]
        Fz_beam[ii] = FReactionHist[end,1]
        M25_beam[ii] = FReactionHist[end,4]
        Mcurv_beam[ii] = FReactionHist[end,5]
        Msweep_beam[ii] = FReactionHist[end,6]
        Ux_beam[ii] = uHist[ii,end-5,end]
        Uy_beam[ii] = uHist[ii,end-4,end]
        Uz_beam[ii] = uHist[ii,end-3,end]
        Θ25_beam[ii] = uHist[ii,end-2,end]
        Θcurv_beam[ii] = uHist[ii,end-1,end]
        Θsweep_beam[ii] = uHist[ii,end,end]

    end #for load

    # for ilevel = 1:length(RPM)
    #     println(Ux_beam[ilevel])
    #     println(Uy_beam[ilevel])
    #     println(Uz_beam[ilevel])
    # end
    # iload idof, itime

    meshx = mesh.x
    meshy = mesh.y
    meshz = mesh.z

    pointvec = zeros(length(points),3)
    for (i,point) in enumerate(points)
        pointvec[i,:] = point[:]
    end

    return meshx, meshy, meshz,pointvec, deformedxyz,dispAnalytical,Ux_beam,uHist,RPM,Uy_beam,Uz_beam
end

meshx, meshy, meshz, pointvec, deformedxyz,dispAnalytical,Ux_beam, uHist,RPM,Uy_beam,Uz_beam = runme(true)
_, _, _, _, deformedxyz_nl,_,Ux_beam_nl, uHist_nl,_,Uy_beam_nl,Uz_beam_nl = runme(false) #nonlinear


# ##############################################
# ####### PLOT
# ##############################################
# deformFact = 1
#
# import PyPlot
# PyPlot.close("all")
# PyPlot.rc("figure", figsize=(4, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# # rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
#
#
# for iload = 1:length(RPM)
#
#     Ux = uHist[iload,1:6:end,end]
#     Uy = uHist[iload,2:6:end,end]
#     Uz = uHist[iload,3:6:end,end]
#
#     Ux_nl = uHist_nl[iload,1:6:end,end]
#     Uy_nl = uHist_nl[iload,2:6:end,end]
#     Uz_nl = uHist_nl[iload,3:6:end,end]
#
#     PyPlot.figure()
#     PyPlot.title("RPM: $(RPM[iload])")
#     PyPlot.plot(meshx,meshy,"k-",label="Undeformed")
#     PyPlot.plot(pointvec[:,1]+deformedxyz[:,1,iload]*deformFact,pointvec[:,2]+deformedxyz[:,2,iload]*deformFact.+0.0015,"-",color = plot_cycle[2],label="GXBeam Linear")
#     PyPlot.plot((meshx+Ux*deformFact),(meshy+Uy*deformFact),"-",color = plot_cycle[1],label="OWENS Linear")
#     PyPlot.plot(pointvec[:,1]+deformedxyz_nl[:,1,iload]*deformFact,pointvec[:,2]+deformedxyz_nl[:,2,iload]*deformFact.+0.0015,"--",color=plot_cycle[2],label="GXBeam Nonlinear")
#     PyPlot.plot((meshx+Ux_nl*deformFact),(meshy+Uy_nl*deformFact),"--",color=plot_cycle[1],label="OWENS Nonlinear")
#     PyPlot.legend()
#     PyPlot.xlabel("x-position (m)")
#     PyPlot.ylabel("y-position (m)")
#     PyPlot.axis("equal")
#     PyPlot.xlim([0,0.5])
#     # PyPlot.savefig("./rotating_$iload.pdf",transparent = true)
# end
#
# PyPlot.rc("figure.subplot", left=.2, bottom=.19, top=0.9, right=.9)
# PyPlot.figure()
# # PyPlot.plot(RPM,dispAnalytical,label="Anaytical (Linear)")
# PyPlot.plot(RPM,deformedxyz[end,2,:].+0.00013,color=plot_cycle[2],label="GXBeam Linear")
# PyPlot.plot(RPM,Uy_beam,color=plot_cycle[1],label="OWENS Linear")
# PyPlot.plot(RPM,deformedxyz_nl[end,2,:].+0.00013,"--",color=plot_cycle[2],label="GXBeam Nonlinear")
# PyPlot.plot(RPM,Uy_beam_nl,"--",color=plot_cycle[1],label="OWENS Nonlinear")
# PyPlot.xlabel("RPM")
# PyPlot.ylabel("Tip Deflection (m)")
# PyPlot.legend()
# # PyPlot.savefig("./tipDeflection.pdf",transparent = true)

###############################################
######## TEST
###############################################

for iload = 1:length(RPM)
    atol = deformedxyz[end,2,iload]*0.02 # 2%
    @test isapprox(Uy_beam[iload],deformedxyz[end,2,iload];atol)

    atol = deformedxyz_nl[end,2,iload]*0.02 # 2%
    @test isapprox(Uy_beam_nl[iload],deformedxyz_nl[end,2,iload];atol)
end
