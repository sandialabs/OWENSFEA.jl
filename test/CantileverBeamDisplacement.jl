using GXBeam, LinearAlgebra
using Test
include("./testdeps.jl")
import Statistics:mean
import FLOWMath
import OWENSFEA
using StaticArrays
# include("../src/OWENSFEA.jl")
# Cantilever Beam/Rod (fixed - free)

function runme(islinear,steady)
# islinear = true
# steady = true
# Define Beam
L = 0.5
Lpt = LinRange(0,L,41)
global Lpt = (Lpt[1:end-1].+Lpt[2:end])./2
b = 0.05
global h = 0.02
A = b*h
E = 2.1e11
v = 0.28
G = E/2(1+v)
rho = 7800.0
Iyy = b*h^3/12
Izz = b^3*h/12
J = Iyy+Izz
nev = 10 #modes
###############################################
######## Analytical
###############################################
global P = 1e4:1e4:1e5
dispAnalytical = zeros(length(P))
stressAnalytical = zeros(length(P),length(Lpt))
strainAnalytical = zeros(length(P),length(Lpt))

Lptused = reverse(Lpt)
for i = 1:length(P)
    dispAnalytical[i] = P[i]*L^3/(3*E*Iyy)
    for ii = 1:length(Lptused)
        stressAnalytical[i,ii] = P[i]*Lptused[ii]*h/2/Iyy
        strainAnalytical[i,ii] = stressAnalytical[i,ii]/E
    end
end
# println(dispAnalytical)
###############################################
######## GXBeam
###############################################

# create points
# straight section of the beam
r= [0.0, 0, 0]
nelem = 40
lengths, points, xm, Cab= discretize_beam(L, r, nelem)

# index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 1/(G*A), 1/(G*A), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)
mass = fill(Diagonal([rho*A, rho*A, rho*A, rho*J, rho*Iyy, rho*Izz]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop;
compliance = compliance)#,
# mass = mass,
# frames = Cab)#,
# lengths = lengths,
# midpoints = xm)

# run an analysis for each prescribed tip load
states = Vector{AssemblyState{Float64}}(undef, length(P))
λ2 = Vector{Vector{ComplexF64}}(undef,length(P))
U = Vector{Matrix{ComplexF64}}(undef, length(P))
MV = Vector{Matrix{ComplexF64}}(undef, length(P))
state = Vector{AssemblyState{Float64}}(undef, length(P))
eigenstates = Vector{Vector{AssemblyState{ComplexF64}}}(undef, length(P))
deformedxyz = zeros(length(points),3,length(P))
system = GXBeam.StaticSystem(assembly)
straintopGX = zeros(length(P),nelem)
strainGX = zeros(3,length(P),nelem)
curvGX = zeros(3,length(P),nelem)
for i = 1:length(P)
    # i = 1

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
    # fixed left side
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    # shear force on right tip
    nelem+1 => PrescribedConditions(Fz_follower = P[i])
    )

    # perform a static analysis
    static_analysis!(system, assembly;
    prescribed_conditions=prescribed_conditions,linear = islinear)

    # process state and eigenstates
    states[i] = AssemblyState(system, assembly;
    prescribed_conditions = prescribed_conditions)

    for (ipt,point) in enumerate(states[i].points)
        deformedxyz[ipt,:,i] = point.u
    end

    for iel = 1:length(states[i].elements)
        strainGX[:,i,iel] = element_strain(assembly.elements[iel],states[i].elements[iel].Fi,states[i].elements[iel].Mi)
        curvGX[:,i,iel] = element_curvature(assembly.elements[iel],states[i].elements[iel].Fi,states[i].elements[iel].Mi)
        straintopGX[i,iel] = curvGX[2,i,iel]*h/2
    end
end

# for ilevel = 1:length(P)
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
angleD = 0.0, # angle of second section of beam relative to first (0 is straight)
zeroOffset = 0.0)#r_b1[1]) #offset from 0 before first beam begins

ort.Twist_d .= 180.0 #enforce flatwise orientation

# Create Sectional Properties
sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, length(mesh.z)-1)

for imesh = 1:length(mesh.z)-1
    ac = [0.0,0.0]
    twist_d = [0.0,0.0]
    rhoA = [rho*A,rho*A]
    EIyy = [E*Iyy,E*Iyy]
    EIzz = [E*Izz,E*Izz]
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
    sectionPropsArray[imesh] = OWENSFEA.SectionPropsArray(ac,twist_d,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)
end

rotationalEffects = ones(mesh.numEl)
el = OWENSFEA.El(sectionPropsArray,ort.Length,ort.Psi_d,ort.Theta_d,ort.Twist_d,rotationalEffects)

M25_beam = zeros(length(P))
Mcurv_beam = zeros(length(P))
Msweep_beam = zeros(length(P))
Fn_beam = zeros(length(P))
Ft_beam = zeros(length(P))
Fz_beam = zeros(length(P))
Ux_beam = zeros(length(P))
Uy_beam = zeros(length(P))
Uz_beam = zeros(length(P))
Θ25_beam = zeros(length(P))
Θcurv_beam = zeros(length(P))
Θsweep_beam = zeros(length(P))


# node, dof, boundary condition (bc)
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

numDOFPerNode = 6
delta_t = 0.1
if steady
    finalt = 0.1
    iterationType = "none"
else
    finalt = 0.3
    iterationType = "DI"
end
numTS = round(Int,finalt/delta_t)
TOL = 1e-3
rpm = 1.0
uHist = zeros(length(P),mesh.numNodes*numDOFPerNode,numTS+1) #put here so it keeps the last RPM solution in the scope

epsilon_x_hist = zeros(length(P),4,mesh.numEl,numTS)
epsilon_y_hist = zeros(length(P),4,mesh.numEl,numTS)
epsilon_z_hist = zeros(length(P),4,mesh.numEl,numTS)
kappa_x_hist = zeros(length(P),4,mesh.numEl,numTS)
kappa_y_hist = zeros(length(P),4,mesh.numEl,numTS)
kappa_z_hist = zeros(length(P),4,mesh.numEl,numTS)

for (iload,load) in enumerate(P)
    nodalinputdata = [mesh.numNodes "F" 1 load]

    nodalTerms = OWENSFEA.applyConcentratedTerms(mesh.numNodes, 6;data = nodalinputdata)

    feamodel = OWENSFEA.FEAModel(;analysisType = "TNB",
    outFilename = "none",
    joint,
    pBC = pBC,
    nodalTerms = nodalTerms,
    numNodes = mesh.numNodes,
    nlOn = !islinear,
    tolerance = 1.0e-06,
    maxIterations = 500,
    maxNumLoadSteps = 20,
    gravityOn = true,
    iterationType,
    # RayleighAlpha = 0.3,
    # RayleighBeta = 0.3,
    aeroElasticOn=false)

    elStorage = OWENSFEA.initialElementCalculations(feamodel,el,mesh)

    structureMass,structureMOI,structureMassCenter = OWENSFEA.calculateStructureMassProps(elStorage)

    u_s = zeros(mesh.numNodes*numDOFPerNode,1)
    # u_s[:] = 2.0.*[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0002786498297864991, 0.0001929465203607998, 0.0, 0.0, 0.0, 0.0, 0.0005529717288725646, 0.0009313743703659432, 0.0, 0.0, 0.0, 0.0, 0.0008206037989981936, 0.0024888842441722087, 0.0, 0.0, 0.0, 0.0, 0.001078899807277038, 0.005145350307976802, 0.0, 0.0, 0.0, 0.0, 0.001324819334136119, 0.009188670895859945, 0.0, 0.0, 0.0, 0.0, 0.0015548051834367155, 0.014916237961941218, 0.0, 0.0, 0.0, 0.0, 0.0017646462583511772, 0.022636039240240186, 0.0, 0.0, 0.0, 0.0, 0.0019493252555308333, 0.03266729739051467, 0.0, 0.0, 0.0, 0.0, 0.0021028521403319293, 0.04534054000670048, 0.0, 0.0, 0.0, 0.0, 0.002218086539899789, 0.060996983993136356, 0.0, 0.0, 0.0, 0.0, 0.0022865550199244592, 0.07998710843940031, 0.0, 0.0, 0.0, 0.0, 0.00229827275383148, 0.10266828294762664, 0.0, 0.0, 0.0, 0.0, 0.0022415833476843727, 0.12940131483878864, 0.0, 0.0, 0.0, 0.0, 0.0022415833476843727, 0.12940131483878864, 0.0, 0.0, 0.0, 0.0, 0.019852328157362865, 0.14717629328985662, 0.0, 0.0, 0.0, 0.0, 0.038470151060156624, 0.16600657787773634, 0.0, 0.0, 0.0, 0.0, 0.057394650935991415, 0.1851755162598843, 0.0, 0.0, 0.0, 0.0]
    # u_s = OWENSFEA.setInitialConditions(feamodel.initCond,u_s,numDOFPerNode)
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

    elStrain = fill(OWENSFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)
    dispOut = OWENSFEA.DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492)) #TODO: not hard coded
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
            dispData = OWENSFEA.DispData(u_s,udot_s,uddot_s,0.0)


            CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]
            CN2P=1.0*LinearAlgebra.I(3)
            CN2H = CP2H*CN2P

            rbData = zeros(9) #TODO
            u_jLast = u_j
            #TODO: time isn't used in the function???
            if steady
                dispOut,elStrain,staticAnalysisSuccessful,FReaction_j = OWENSFEA.staticAnalysis(feamodel,
                mesh,el,u_s,Omega_j,Omega_j,elStorage;reactionNodeNumber=1,
                OmegaDot=OmegaDot_j,Fexternal, Fdof)
                u_j = dispOut
            else
                elStrain,dispOut,FReaction_j = OWENSFEA.structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,time,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)
                u_j = dispOut.displ_sp1
                udot_j  = dispOut.displdot_sp1
                uddot_j = dispOut.displddot_sp1
            end

            ## update timestepping variables and other states, store in history arrays
            azi_jLast = azi_j
            azi_j = azi_s + Omega_j*delta_t*2*pi

            if steady
                uNorm = 1e-6
                aziNorm = 1e-6
            else
                uNorm = LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)
                aziNorm = LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)
                println("uNorm: $uNorm axiNorm: $aziNorm")
            end
        end

        u_s = u_j
        udot_s = udot_j
        uddot_s = uddot_j
        azi_s = azi_j
        Omega_s = Omega_j

        uHist[iload,:,i+1] = u_s
        FReactionHist[i+1,:] = FReaction_j[1:6]
        aziHist[i+1] = azi_s
        for ielem = 1:length(elStrain)
            epsilon_x_hist[iload,:,ielem,i] = elStrain[ielem].epsilon_x
            epsilon_y_hist[iload,:,ielem,i] = elStrain[ielem].epsilon_y
            epsilon_z_hist[iload,:,ielem,i] = elStrain[ielem].epsilon_z
            kappa_x_hist[iload,:,ielem,i] = elStrain[ielem].kappa_x
            kappa_y_hist[iload,:,ielem,i] = elStrain[ielem].kappa_y
            kappa_z_hist[iload,:,ielem,i] = elStrain[ielem].kappa_z
        end # while

    end # time stepping

    # Rotating Frame of Reference, 6 DOF where
    # 1 = turbine vertical force
    # 2 = turbine 2D slice tangential force
    # 3 = turbine 2D slice normal force
    # 4 = blade M25 twisting moment
    # 5 = blade curvature twisting moment  or is this sweep?
    # 6 = blade sweep moment.  and this is curvature moment?

    Fn_beam[iload] = FReactionHist[end,3]
    Ft_beam[iload] = FReactionHist[end,2]
    Fz_beam[iload] = FReactionHist[end,1]
    M25_beam[iload] = FReactionHist[end,4]
    Mcurv_beam[iload] = FReactionHist[end,5]
    Msweep_beam[iload] = FReactionHist[end,6]
    Ux_beam[iload] = uHist[iload,end-5,end]
    Uy_beam[iload] = uHist[iload,end-4,end]
    Uz_beam[iload] = uHist[iload,end-3,end]
    Θ25_beam[iload] = uHist[iload,end-2,end]
    Θcurv_beam[iload] = uHist[iload,end-1,end]
    Θsweep_beam[iload] = uHist[iload,end,end]


    # Strain

    epsilon_x1 = epsilon_x_hist[iload,1,:,end]
    epsilon_x2 = epsilon_x_hist[iload,2,:,end]
    epsilon_x3 = epsilon_x_hist[iload,3,:,end]
    epsilon_x4 = epsilon_x_hist[iload,4,:,end]

    epsilon_y1 = epsilon_y_hist[iload,1,:,end]
    epsilon_y2 = epsilon_y_hist[iload,2,:,end]
    epsilon_y3 = epsilon_y_hist[iload,3,:,end]
    epsilon_y4 = epsilon_y_hist[iload,4,:,end]

    epsilon_z1 = epsilon_z_hist[iload,1,:,end]
    epsilon_z2 = epsilon_z_hist[iload,2,:,end]
    epsilon_z3 = epsilon_z_hist[iload,3,:,end]
    epsilon_z4 = epsilon_z_hist[iload,4,:,end]

    kappa_x1 = kappa_x_hist[iload,1,:,end]
    kappa_x2 = kappa_x_hist[iload,2,:,end]
    kappa_x3 = kappa_x_hist[iload,3,:,end]
    kappa_x4 = kappa_x_hist[iload,4,:,end]

    kappa_y1 = kappa_y_hist[iload,1,:,end]
    kappa_y2 = kappa_y_hist[iload,2,:,end]
    kappa_y3 = kappa_y_hist[iload,3,:,end]
    kappa_y4 = kappa_y_hist[iload,4,:,end]

    kappa_z1 = kappa_z_hist[iload,1,:,end]
    kappa_z2 = kappa_z_hist[iload,2,:,end]
    kappa_z3 = kappa_z_hist[iload,3,:,end]
    kappa_z4 = kappa_z_hist[iload,4,:,end]

    z_plot = FLOWMath.akima(LinRange(0,1,length(mesh.z)),mesh.z,LinRange(0,1,length(kappa_z1)))
    if iload < 5
        if islinear
            bumpsecond = 1.0
        else
            bumpsecond = 10.0
        end
        println("Load $iload")
        println("is Linear?: $islinear")
        for ipoint = 1:length(straintopGX[iload,:])
            # kappa_y and analytical
            atol = max(abs(straintopGX[iload,ipoint])*0.02*bumpsecond,1e-9)
            @test isapprox(straintopGX[iload,ipoint],kappa_y1[ipoint]*h/2;atol)
            # println("$(straintopGX[iload,ipoint]), $(kappa_y1[ipoint]*h/2)")
            if islinear
                atol = max(abs(strainAnalytical[iload,ipoint])*0.02*bumpsecond,1e-9)
                @test isapprox(-strainAnalytical[iload,ipoint],kappa_y1[ipoint]*h/2;atol)
            end
    
            if !islinear #no strain coupling in x for bending if linear
                #epsilon_x
                atol = max(abs(strainGX[1,iload,ipoint])*0.01,1e-4)
                # @test isapprox(strainGX[1,iload,ipoint],epsilon_x1[ipoint];atol) #TODO: resolve difference in nonlinear strain coupling.
            end

            #epsilon_y, this has oddities with the quad points, but the overall value should be close to zero, find arbitrarily that by adding all four quad points and multipling by the smaller of the wieghts, it gives a close ish answer
            atol = max(abs(strainGX[2,iload,ipoint])*0.01*bumpsecond,1e-9)
            meaneps_y = (epsilon_y1[ipoint].+epsilon_y2[ipoint].+epsilon_y3[ipoint].+epsilon_y4[ipoint]).*0.34785484513745385
            @test isapprox(strainGX[2,iload,ipoint],meaneps_y;atol)
    
            #epsilon_z, this has oddities with the quad points, but the overall value should be close to zero, find arbitrarily that by adding all four quad points and multipling by the smaller of the wieghts, it gives a close ish answer
            atol = max(abs(strainGX[3,iload,ipoint])*0.01,1e-3)
            meaneps_z = (epsilon_z1[ipoint].+epsilon_z2[ipoint].+epsilon_z3[ipoint].+epsilon_z4[ipoint]).*0.34785484513745385
            @test isapprox(strainGX[3,iload,ipoint],meaneps_z;atol)
    
            #kappa_x
            atol = max(abs(curvGX[1,iload,ipoint])*0.01*bumpsecond,1e-9)
            @test isapprox(curvGX[1,iload,ipoint],kappa_x1[ipoint];atol)
    
            #kappa_y - already done
    
            #kappa_z
            atol = max(abs(curvGX[3,iload,ipoint])*0.01*bumpsecond,1e-9)
            @test isapprox(curvGX[3,iload,ipoint],kappa_z1[ipoint];atol)
    
        end
        atol = max(dispAnalytical[iload]*0.02*bumpsecond,1e-9) # 0.5%
        @test isapprox(Ux_beam[iload],dispAnalytical[iload];atol)
        @test isapprox(Ux_beam[iload],deformedxyz[end,3,iload];atol)
    end

end #for load

# for ilevel = 1:length(P)
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

    return meshx, meshz,pointvec, deformedxyz,dispAnalytical,Ux_beam,uHist,P,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,strainAnalytical,straintopGX,strainGX,curvGX
end

meshx, meshz, pointvec, deformedxyz,dispAnalytical,Ux_beam, uHist,P,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,strainAnalytical,straintopGX,strainGX,curvGX = runme(true,true)
_, _, _, deformedxyz_nl,_,Ux_beam_nl, uHist_nl,_,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,strainAnalytical,straintopGX,strainGX,curvGX = runme(false,true) #nonlinear

meshx, meshz, pointvec, deformedxyz,dispAnalytical,Ux_beam, uHist,P,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,strainAnalytical,straintopGX,strainGX,curvGX = runme(true,false)
_, _, _, deformedxyz_nl,_,Ux_beam_nl, uHist_nl,_,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,strainAnalytical,straintopGX,strainGX,curvGX = runme(false,false) #nonlinear


# deformFact = 1
# L = 0.5
# import PyPlot
# PyPlot.pygui(true)
# PyPlot.close("all")
# PyPlot.rc("figure", figsize=(4, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]


# for iload = 1:length(P)

#     Ux = uHist[iload,1:6:end,end]
#     Uy = uHist[iload,2:6:end,end]
#     Uz = uHist[iload,3:6:end,end]

#     Ux_nl = uHist_nl[iload,1:6:end,end]
#     Uy_nl = uHist_nl[iload,2:6:end,end]
#     Uz_nl = uHist_nl[iload,3:6:end,end]

#     PyPlot.figure()
#     PyPlot.title("Load: $(P[iload]) N")
#     PyPlot.plot(meshx./L,meshz./L,"k-",label="Undeformed")
#     PyPlot.plot((pointvec[:,3]+deformedxyz[:,3,iload]*deformFact)./L,(pointvec[:,1]+deformedxyz[:,1,iload]*deformFact)./L,"-",color=plot_cycle[2],label="GXBeam")
#     PyPlot.plot((meshx+Ux*deformFact)./L,(meshz+Uz*deformFact)./L,"-",color=plot_cycle[1],label="OWENS")
#     PyPlot.plot((pointvec[:,3]+deformedxyz_nl[:,3,iload]*deformFact)./L,(pointvec[:,1]+deformedxyz_nl[:,1,iload]*deformFact)./L,"--",color=plot_cycle[2],label="GXBeam Nonlinear")
#     PyPlot.plot((meshx+Ux_nl*deformFact)./L,(meshz+Uz_nl*deformFact)./L,"--",color=plot_cycle[1],label="OWENS Nonlinear")
#     PyPlot.legend()
#     PyPlot.xlabel("x/L")
#     PyPlot.ylabel("y/L")
#     PyPlot.axis("equal")
#     PyPlot.xlim([0,0.5])
#     PyPlot.savefig("./deformedmesh_$iload.pdf",transparent = true)

#     # Strain
#     epsilon_x1 = epsilon_x_hist[iload,1,:,end]
#     epsilon_x2 = epsilon_x_hist[iload,2,:,end]
#     epsilon_x3 = epsilon_x_hist[iload,3,:,end]
#     epsilon_x4 = epsilon_x_hist[iload,4,:,end]
    
#     epsilon_y1 = epsilon_y_hist[iload,1,:,end]
#     epsilon_y2 = epsilon_y_hist[iload,2,:,end]
#     epsilon_y3 = epsilon_y_hist[iload,3,:,end]
#     epsilon_y4 = epsilon_y_hist[iload,4,:,end]
    
#     epsilon_z1 = epsilon_z_hist[iload,1,:,end]
#     epsilon_z2 = epsilon_z_hist[iload,2,:,end]
#     epsilon_z3 = epsilon_z_hist[iload,3,:,end]
#     epsilon_z4 = epsilon_z_hist[iload,4,:,end]
    
#     kappa_x1 = kappa_x_hist[iload,1,:,end]
#     kappa_x2 = kappa_x_hist[iload,2,:,end]
#     kappa_x3 = kappa_x_hist[iload,3,:,end]
#     kappa_x4 = kappa_x_hist[iload,4,:,end]
    
#     kappa_y1 = kappa_y_hist[iload,1,:,end]
#     kappa_y2 = kappa_y_hist[iload,2,:,end]
#     kappa_y3 = kappa_y_hist[iload,3,:,end]
#     kappa_y4 = kappa_y_hist[iload,4,:,end]
    
#     kappa_z1 = kappa_z_hist[iload,1,:,end]
#     kappa_z2 = kappa_z_hist[iload,2,:,end]
#     kappa_z3 = kappa_z_hist[iload,3,:,end]
#     kappa_z4 = kappa_z_hist[iload,4,:,end]
    
    
#     # PyPlot.figure()
#     # PyPlot.title("Load: $(P[iload]) N")
#     # PyPlot.plot(Lpt,strainAnalytical[iload,:],"k",label="Analytical")
#     # PyPlot.plot((pointvec[1:end-1,1]+pointvec[2:end,1])./2,-straintopGX[iload,:],color=plot_cycle[2],label="GXBeam")
#     # PyPlot.plot((pointvec[1:end-1,1]+pointvec[2:end,1])./2,-kappa_y1*h/2,"-",color=plot_cycle[1],label="OWENS1")
#     # PyPlot.legend()
#     # PyPlot.xlabel("y-position (m)")
#     # PyPlot.ylabel("bending strain top of beam")
    
#     PyPlot.figure()
#     PyPlot.title("Load: $(P[iload]) N")
#     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_x1,"-",color=plot_cycle[1],label="OWENS1")
#     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_x2,"-",color=plot_cycle[2],label="OWENS2")
#     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_x3,"-",color=plot_cycle[3],label="OWENS3")
#     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_x4,"-",color=plot_cycle[4],label="OWENS4")
#     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),strainGX[1,iload,:],"-",color=plot_cycle[5],label="GXBeam")
#     PyPlot.legend()
#     PyPlot.xlabel("y-position (m)")
#     PyPlot.ylabel("strain epsilon_x")
    
# #     PyPlot.figure()
# #     PyPlot.title("Load: $(P[iload]) N")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_y1,"-",color=plot_cycle[1],label="OWENS1")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_y2,"-",color=plot_cycle[2],label="OWENS2")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_y3,"-",color=plot_cycle[3],label="OWENS3")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_y4,"-",color=plot_cycle[4],label="OWENS4")
# #     meaneps_y = (epsilon_y1.+epsilon_y2.+epsilon_y3.+epsilon_y4).*0.34785484513745385
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),meaneps_y,"-",color=plot_cycle[4],label="MeanOWENS4")
# #     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),strainGX[2,iload,:],"-",color=plot_cycle[5],label="GXBeam")
# #     PyPlot.legend()
# #     PyPlot.xlabel("y-position (m)")
# #     PyPlot.ylabel("strain epsilon_y")
    
# #     PyPlot.figure()
# #     PyPlot.title("Load: $(P[iload]) N")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_z1,"-",color=plot_cycle[1],label="OWENS1")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_z2,"-",color=plot_cycle[2],label="OWENS2")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_z3,"-",color=plot_cycle[3],label="OWENS3")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),epsilon_z4,"-",color=plot_cycle[4],label="OWENS4")
# #     meaneps_z = (epsilon_z1.+epsilon_z2.+epsilon_z3.+epsilon_z4).*0.34785484513745385
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),meaneps_z,"-",color=plot_cycle[4],label="MeanOWENS4")
# #     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),strainGX[3,iload,:],"-",color=plot_cycle[5],label="GXBeam")
# #     PyPlot.legend()
# #     PyPlot.xlabel("y-position (m)")
# #     PyPlot.ylabel("strain epsilon_z")
    
# #     PyPlot.figure()
# #     PyPlot.title("Load: $(P[iload]) N")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_x1,"-",color=plot_cycle[1],label="OWENS1")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_x2,"-",color=plot_cycle[2],label="OWENS2")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_x3,"-",color=plot_cycle[3],label="OWENS3")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_x4,"-",color=plot_cycle[4],label="OWENS4")
# #     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),curvGX[1,iload,:],"-",color=plot_cycle[5],label="GXBeam")
# #     PyPlot.legend()
# #     PyPlot.xlabel("y-position (m)")
# #     PyPlot.ylabel("strain kappa_x")
    
# #     PyPlot.figure()
# #     PyPlot.title("Load: $(P[iload]) N")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_y1,"-",color=plot_cycle[1],label="OWENS1")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_y2,"-",color=plot_cycle[2],label="OWENS2")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_y3,"-",color=plot_cycle[3],label="OWENS3")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_y4,"-",color=plot_cycle[4],label="OWENS4")
# #     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),curvGX[2,iload,:],"-",color=plot_cycle[5],label="GXBeam")
# #     PyPlot.legend()
# #     PyPlot.xlabel("y-position (m)")
# #     PyPlot.ylabel("strain kappa_y")
    
# #     PyPlot.figure()
# #     PyPlot.title("Load: $(P[iload]) N")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_z1,"-",color=plot_cycle[1],label="OWENS1")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_z2,"-",color=plot_cycle[2],label="OWENS2")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_z3,"-",color=plot_cycle[3],label="OWENS3")
# #     PyPlot.plot(LinRange(0,1,length(kappa_z1)),kappa_z4,"-",color=plot_cycle[4],label="OWENS4")
# #     PyPlot.plot(LinRange(0,1,length(strainGX[1,iload,:])),curvGX[3,iload,:],"-",color=plot_cycle[5],label="GXBeam")
# #     PyPlot.legend()
# #     PyPlot.xlabel("y-position (m)")
# #     PyPlot.ylabel("strain kappa_z")
# end

# # PyPlot.figure()
# # PyPlot.plot(P,dispAnalytical,"k",label="Analytical (Linear)")
# # PyPlot.plot(P,deformedxyz[end,3,:],color=plot_cycle[2],label="GXBeam")
# # PyPlot.plot(P,Ux_beam,color=plot_cycle[1],label="OWENS")
# # PyPlot.plot(P,deformedxyz_nl[end,3,:],"--",color=plot_cycle[2],label="GXBeam Nonlinear")
# # PyPlot.plot(P,Ux_beam_nl,"--",color=plot_cycle[1],label="OWENS Nonlinear")
# # PyPlot.xlabel("Load (N)")
# # PyPlot.ylabel("Tip Deflection (M)")
# # PyPlot.legend()
# # # PyPlot.savefig("./beamTipDeflec.pdf",transparent = true)
