using GXBeam, LinearAlgebra
import GyricFEA
# include("../src/GyricFEA.jl")
include("./testdeps.jl")
import Statistics
using Test
# Cantilever Beam/Rod (fixed - free)

# Define Beam
L = 0.5
b = 0.05
h = 0.02#02
A = b*h
E = 2.1e11
v = 0.28
G = E/2(1+v)
rho = 7800.0
Iyy = b*h^3/12
Izz = b^3*h/12
J = Iyy+Izz
RPM = 0:1000:6000
nev = 12 #modes
sweepD = 45.0

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
# points = [[points1[i][3],points1[i][2],points1[i][1]] for i = 1:length(points1)]

start = 1:nelem_b1 + nelem_b2
stop = 2:nelem_b1 + nelem_b2 + 1
lengths = vcat(lengths_b1, lengths_b2)
xm = vcat(xm_b1, xm_b2)
# xm = [[xm1[i][3],xm1[i][2],xm1[i][1]] for i = 1:length(xm1)]
Cab = vcat(Cab_b1, Cab_b2)

# compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 1/(E*A/2.6*5/6), 1/(E*A/2.6*5/6), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)
mass = fill(Diagonal([rho*A, rho*A, rho*A, rho*J, rho*Iyy, rho*Izz]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop;
    compliance = compliance,
    mass = mass,
    frames = Cab,
    lengths = lengths,
    midpoints = xm)

system = System(assembly,false)
# create dictionary of prescribed conditions
prescribed_conditions = Dict(
# fixed left side
1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
# shear force on right tip
# nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
)

freqGXBeam = zeros(length(RPM),Int(nev/2))
for j = 1:length(RPM)
# j = 2
    w0 = [0, 0,RPM[j]*(2*pi)/60]

    # eigenvalues and (right) eigenvectors
    global system
    system, λ2, V, converged = eigenvalue_analysis!(system, assembly;
    prescribed_conditions = prescribed_conditions,
    angular_velocity = w0,
    linear = false,
    reset_state = true,
    find_steady_state = true,
    nev = nev)

    # corresponding left eigenvectors
    U = left_eigenvectors(system, λ2, V)

    # post-multiply mass matrix with right eigenvector matrix
    # (we use this later for correlating eigenvalues)
    MV = system.M * V

    # process state and eigenstates
    state = AssemblyState(system, assembly;
    prescribed_conditions = prescribed_conditions)
    eigenstates = [AssemblyState(system, assembly, V[:,k];
        prescribed_conditions = prescribed_conditions) for k = 1:nev]

    frequencyNative = [imag(λ2[k])/(2*pi) for k = 1:2:nev]

    # set previous left eigenvector matrix
    U_p = copy(U)

    # construct correlation matrix
    C = U_p*MV

    # correlate eigenmodes
    perm, corruption = correlate_eigenmodes(C)

    # re-arrange eigenvalues and eigenvectors
    λ2 = λ2[perm]
    # U = U[perm,:]
    # MV = MV[:,perm]
    # eigenstates = eigenstates[perm]
    #
    # # update previous eigenvector matrix
    # U_p .= U
    #
    # # update previous eigenvector matrix
    # U_p .= U[1]

    freqGXBeam[j,:] = [imag(λ2[k])/(2*pi) for k = 1:2:nev]
end
# println(freqGXBeam)

#These should match the FEA case from https://autofem.com/examples/determining_natural_frequencie.html
# with frequencies of 67, 418, 1157

#NOTE: with the full stiffness and mass matrices the GX modes matching the 2D analytical in plane modes are 1,3,5.

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

# node, dof, boundary condition (bc)
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]
# mesh.numNodes 1 0
# mesh.numNodes 2 0
# mesh.numNodes 3 0
# mesh.numNodes 4 0
# mesh.numNodes 5 0
# mesh.numNodes 6 0]

feamodel = GyricFEA.FEAModel(;analysisType = "M",
outFilename = "none",
joint = joint,
gravityOn = false,
platformTurbineConnectionNodeNumber = 1,
pBC = pBC,
nlOn = true,
spinUpOn = true,
numNodes = mesh.numNodes)

freqOWENS = zeros(length(RPM),Int(nev/2))
for j = 1:length(RPM)
    # println(j)
    Omega = RPM[j]/60 #*(2*pi)
    OmegaStart = RPM[j]/60 #*(2*pi)

    freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=GyricFEA.modal(feamodel,mesh,el;Omega,OmegaStart)
    # OWENS Frequencies that correspond to the GX beam are every other, and then 1,3,5 of the every other sets corresponds to the analytical
    freqOWENS[j,:] = freq[1:2:nev]
end

###############################################
######## TEST
###############################################
# println(freqOWENS)
# println(freqGXBeam)

for i = 1:length(RPM)
    for j = 1:Int(nev/2)
        # println(i)
        if i>=5
            atol = freqGXBeam[i,j]*0.05
        else
            atol = freqGXBeam[i,j]*0.03
        end
        @test isapprox(freqGXBeam[i,j], freqOWENS[i,j];atol)
    end
end
#
# ###############################################
# ######## PLOT
# ###############################################
# import PyPlot
# PyPlot.ion()
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
# PyPlot.figure()
# for i=1:3:18
#     linex=[RPM[1], RPM[end]+10]
#     liney=[RPM[1], RPM[end]+10].*i./60.0
#     PyPlot.plot(linex,liney,"--k",linewidth=0.5)
#     PyPlot.annotate("$i P",xy=(0.95*linex[2],liney[2]+.05+(i-1)*.01))
# end
# for i = 1:Int(nev/2)
#     PyPlot.plot(RPM,freqGXBeam[:,i],color=plot_cycle[2])
#     PyPlot.plot(RPM,freqOWENS[:,i],color=plot_cycle[1])
# end
# PyPlot.plot(0,0,color=plot_cycle[2],label="GXBeam")
# PyPlot.plot(0,0,color=plot_cycle[1],label="OWENS")
# PyPlot.xlabel("RPM")
# PyPlot.ylabel("Frequency (Hz)")
# PyPlot.legend(loc=(0.05,0.68))
# PyPlot.savefig("rotating_modal.pdf",transparent = true)
