using GXBeam, LinearAlgebra
import GyricFEA
using Test
include("./testdeps.jl")
# include("../src/GyricFEA.jl")
timestart = time()

L = 60 # m

# create points
nelem = 10
x = range(0, L, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# stiffness matrix for each beam element
# stiffness = fill(
# [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
# 1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
# 6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
# -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
# -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
# -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
# nelem)

stiffness = fill(Diagonal([2.389e9,4.334e8,2.743e7,2.167e7,1.970e7,4.406e8]), nelem)

# mass matrix for each beam element
# mass = fill(
# [258.053      0.0        0.0      0.0      7.07839  -71.6871
# 0.0      258.053      0.0     -7.07839  0.0        0.0
# 0.0        0.0      258.053   71.6871   0.0        0.0
# 0.0       -7.07839   71.6871  48.59     0.0        0.0
# 7.07839    0.0        0.0      0.0      2.172      0.0
# -71.6871     0.0        0.0      0.0      0.0       46.418],
# nelem)

mass = fill(Diagonal([258.053,258.053,258.053,48.59,2.172,46.418]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

# simulation time
dt = 0.001
delta_t = dt
t = 0:dt:2.0
nstep = length(t)

# prescribed conditions
prescribed_conditions = (t) -> begin
    Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # force on right side
        nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t))
        )
end

system, history, converged = time_domain_analysis(assembly, t; prescribed_conditions=prescribed_conditions)

elapsedGX = time() - timestart

function runOWENS()
    ###############################################
    ######## OWENS
    ###############################################
    islinear = false
    # Create Mesh
    mesh, ort, joint = mesh_beam(;L1 = L/2, #first section of beam length
    L2 = L/2, #second section of beam length
    Nelem1 = 5,
    Nelem2 = 5,
    angleD = 0.0, # angle of second section of beam relative to first (0 is straight)
    zeroOffset = 0.0,
    vertical = true)#r_b1[1]) #offset from 0 before first beam begins

    # Create Sectional Properties
    sectionPropsArray = Array{GyricFEA.SectionPropsArray, 1}(undef, length(mesh.z)-1)

    for ii = 1:length(mesh.z)-1
        ac = [0.0,0.0]
        twist_d = [0.0,0.0]
        rhoA = [mass[1][1,1],mass[1][1,1]]
        EIyy = [stiffness[1][5,5],stiffness[1][5,5]]
        EIzz = [stiffness[1][6,6],stiffness[1][6,6]]
        GJ = [stiffness[1][4,4],stiffness[1][4,4]]
        EA = [stiffness[1][1,1],stiffness[1][1,1]]
        rhoIyy = [mass[1][5,5],mass[1][5,5]]
        rhoIzz = [mass[1][6,6],mass[1][6,6]]
        rhoJ = [mass[1][4,4],mass[1][4,4]]
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

    # node, dof, boundary condition (bc)
    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0]

    numDOFPerNode = 6
    numTS = length(t)
    TOL = 1e-3
    rpm = 1.0
    uHist = zeros(mesh.numNodes*numDOFPerNode,numTS+1) #put here so it keeps the last RPM solution in the scope

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
    iterationType="NR",
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

    rpm = 0.0
    Omega_s = rpm/60.0#*2*pi
    Omega_j = Omega_s
    OmegaDot_j = 0.0
    Fexternal = []
    Fdof = []

    elStrain = fill(GyricFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)
    dispOut = GyricFEA.DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492)) #TODO: not hard coded
    FReactionHist = zeros(numTS+1,6)
    aziHist = zeros(numTS+1)
    eps_xx_0_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_z_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_z_hist = zeros(4,mesh.numEl,numTS)

    FReactionsm1 = zeros(6)
    FReactionHist[1,:] = FReactionsm1
    FReaction_j = FReactionsm1

    # Run unsteady analysis
    for (i,time) in enumerate(t)
        # println("Time: $time")
        # println("$i of $(length(t))")
        iwhile = 0
        uNorm = 1e5
        aziNorm = 1e5
        while (uNorm > TOL || aziNorm > TOL) && iwhile < 10 #TODO: add TOL to feamodel
            # println("While Loop Iteration: $iwhile")

            iwhile += 1
            dispData = GyricFEA.DispData(u_s,udot_s,uddot_s,0.0)


            CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]
            CN2P=1.0*LinearAlgebra.I(3)
            CN2H = CP2H*CN2P

            rbData = zeros(9) #TODO

            Fexternal = -1e5*sin(20*time)
            Fdof = mesh.numNodes*6-4 #x-direction when beam is vertical

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
            # println("uNorm: $uNorm axiNorm: $aziNorm")
        end

        u_s = u_j
        udot_s = udot_j
        uddot_s = uddot_j
        azi_s = azi_j
        Omega_s = Omega_j

        uHist[:,i+1] = u_s
        FReactionHist[i+1,:] = FReaction_j
        aziHist[i+1] = azi_s
        for ii = 1:length(elStrain)
            eps_xx_0_hist[:,ii,i] = elStrain[ii].eps_xx_0
            eps_xx_z_hist[:,ii,i] = elStrain[ii].eps_xx_z
            eps_xx_y_hist[:,ii,i] = elStrain[ii].eps_xx_y
            gam_xz_0_hist[:,ii,i] = elStrain[ii].gam_xz_0
            gam_xz_y_hist[:,ii,i] = elStrain[ii].gam_xz_y
            gam_xy_0_hist[:,ii,i] = elStrain[ii].gam_xy_0
            gam_xy_z_hist[:,ii,i] = elStrain[ii].gam_xy_z
        end # while

    end # time stepping

    # Rotating Frame of Reference, 6 DOF where
    # 1 = turbine vertical force
    # 2 = turbine 2D slice tangential force
    # 3 = turbine 2D slice normal force
    # 4 = blade M25 twisting moment
    # 5 = blade curvature twisting moment  or is this sweep?
    # 6 = blade sweep moment.  and this is curvature moment?

    Fn_beam = FReactionHist[:,3]
    Ft_beam = FReactionHist[:,2]
    Fz_beam = FReactionHist[:,1]
    M25_beam = FReactionHist[:,4]
    Mcurv_beam = FReactionHist[:,5]
    Msweep_beam = FReactionHist[:,6]
    Ux_beam = uHist[end-5,:]
    Uy_beam = uHist[end-4,:]
    Uz_beam = uHist[end-3,:]
    Θ25_beam = uHist[end-2,:]
    Θcurv_beam = uHist[end-1,:]
    Θsweep_beam = uHist[end,:]
    return Fn_beam,Ft_beam,Fz_beam,M25_beam,Mcurv_beam,Msweep_beam,Ux_beam,Uy_beam,Uz_beam,Θ25_beam,Θcurv_beam,Θsweep_beam
end

timestart = time()
Fn_beam,Ft_beam,Fz_beam,M25_beam,Mcurv_beam,Msweep_beam,Ux_beam,Uy_beam,Uz_beam,Θ25_beam,Θcurv_beam,Θsweep_beam = runOWENS()
elapsedOW = time() - timestart

########################################
# Plot
########################################

println("Time GXBeam: $elapsedGX")
println("Time OWENS: $elapsedOW")
println("GXBeam is $(elapsedOW/elapsedGX) x faster than OWENS")
# import PyPlot
# PyPlot.ion()
# PyPlot.rc("figure", figsize=(4, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
point = vcat(fill(nelem+1, 6), fill(1, 6))
field = [:u, :u, :u, :theta, :theta, :theta, :F, :F, :F, :M, :M, :M]
direction = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
ylabel = ["\$u_x\$ (\$m\$)", "\$u_y\$ (\$m\$)", "\$u_z\$ (\$m\$)",
"Rodriguez Parameter \$\\theta_x\$ (degree)",
"Rodriguez Parameter \$\\theta_y\$ (degree)",
"Rodriguez Parameter \$\\theta_z\$ (degree)",
"\$F_x\$ (\$N\$)", "\$F_y\$ (\$N\$)", "\$F_z\$ (\$N\$)",
"\$M_x\$ (\$Nm\$)", "\$M_y\$ (\$Nm\$)", "\$M_z\$ (\$N\$)"]

# PyPlot.close("all")

for i = 7:12
    # i = 7
    # PyPlot.figure()
    y = [getproperty(state.points[point[i]], field[i])[direction[i]] for state in history]

    if field[i] == :theta
        # convert to Rodriguez parameter
        y .= 4*atan.(y./4)
        # convert to degrees
        y .= rad2deg.(y)
    end

    if field[i] == :F || field[i] == :M
        y .= -y
    end

    # PyPlot.plot(t, y, color = plot_cycle[2], label = "GXBeam")
    if i == 7
        # PyPlot.plot(t,-Fn_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# N")
        # PyPlot.ylim([-7e6,7e6])
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    elseif i == 8
        # PyPlot.plot(t,Fz_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# Z")
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    elseif i == 9
        # PyPlot.plot(t,Ft_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# T")
        # PyPlot.ylim([-3e6,3e6])
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    elseif i == 10
        # PyPlot.plot(t,Mcurv_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# Mcurv")
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    elseif i == 11
        # PyPlot.plot(t,M25_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# M25")
        # PyPlot.ylim([-4e6,4e6])
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    elseif i == 12
        # PyPlot.plot(t,Msweep_beam[1:end-1], color = plot_cycle[1], label = "OWENS")# Msweep")
        myerror = sum(abs.(-Fn_beam[1:end-1]-y))./sum(abs.(y))
    end
    # PyPlot.xlim([0, 2.0])
    # PyPlot.xticks(collect(0:0.5:2.0))
    # PyPlot.xlabel("Time (s)")
    # PyPlot.ylabel(ylabel[i])
    # PyPlot.legend()
    # PyPlot.savefig("./Unsteady$(ylabel[i]).pdf",transparent = true)

    if i == 7 || i == 9 || i == 11 # the nonzero fx fz my
        @test myerror < 3.3
    end

end
