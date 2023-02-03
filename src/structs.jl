mutable struct FEAModel
    analysisType
    initCond
    aeroElasticOn
    guessFreq
    airDensity
    gravityOn
    nlOn
    spinUpOn
    outFilename
    RayleighAlpha
    RayleighBeta
    elementOrder
    joint
    platformTurbineConnectionNodeNumber
    jointTransform
    reducedDOFList
    nlParams
    BC
    nodalTerms
    numModes
    return_all_reaction_forces
end

# this way you can use defaults and pass in what is different, and it's mapped
# by keyword so it doesn't have to be in order.

"""
    FEAModel(;analysisType = "TNB",
        initCond = [],
        aeroElasticOn = false,
        guessFreq = 0.0,
        airDensity=1.2041,
        gravityOn = true,
        nlOn = false,
        spinUpOn = false,
        outFilename = "none",
        RayleighAlpha = 0.0,
        RayleighBeta = 0.0,
        elementOrder = 1,
        joint = [0,0],
        platformTurbineConnectionNodeNumber = 1,
        jointTransform = 0.0,
        reducedDOFList = zeros(Int,2),
        numDOFPerNode = 6,
        numNodes = 0,
        numModes = 20,
        nlParams = 0,
        pBC = 0,
        nodalTerms = 0.0,
        iterationType = "NR",
        adaptiveLoadSteppingFlag = true,
        tolerance = 1.0000e-06,
        maxIterations = 50,
        maxNumLoadSteps = 20,
        minLoadStepDelta = 0.0500,
        minLoadStep = 0.0500,
        prescribedLoadStep = 0.0,
        return_all_reaction_forces = false)

Model inputs for FEA analysis, struct

# Inputs
* `analysisType::string`: Newmark Beta time stepping "TNB", Dean time stepping "TD", modal "M"
* `initCond::Array{<:float}`: Initial conditions Nx3 matrix consisting of nodeNumber, local DOF (1-6), and displacement value
* `aeroElasticOn::Bool`: Include simplified flutter calculataions in the timoshenko element matrices
* `guessFreq::float`: aeroelastic starting guess, only used if aeroElasticOn
* `airDensity::float`: working fluid density
* `gravityOn::Bool orArray{<:float}`: vector of 3 or flag to include distributed gravity acceleration (9.81m/s) in the negative z-direction
* `nlOn::Bool`: flag for solver to calculate deflection induced stiffness changes and associated convergance to the coupled solution
* `spinUpOn::Bool`: flag to perform a static analysis (warm start) prior to performing modal analysis
* `outFilename::string`: /path/to/desired/output/filename if it doesn't exist already it is created, if exists, is overwritten
* `RayleighAlpha::float`: Rayleigh alpha damping used in timoshenko beam damping matrix
* `RayleighBeta::float`: Rayleigh beta damping used in timoshenko beam damping matrix
* `elementOrder::int`: order of element: 1 linear, 2 quadratic
* `joint::Array{<:float}`: jointNumber masterNode slaveNode jointType (0 weld/fixed, 1 pinned, 2 hinge along local "2", 3 hinge about local "1", 4 hinge along "3", 5 rigid bar constraint)  jointMass 0.0 jointPsi jointTheta
* `platformTurbineConnectionNodeNumber::int`: node at which reaction forces are calculated
* `jointTransform`: not used as an input, is calculated, local transform between dependent and active DOFs for nodes associated with a particular joint
* `reducedDOFList::Array{<:int}`: not used as an input, is calculated, map of original DOF numbering to reduced DOF numbering
* `numDOFPerNode::int`: number of degrees of freedom per node
* `numNodes::int`: total number of nodes in the mesh
* `numModes::int`: number of modes to calculate
* `nlParams::NlParams`: optional there in case the Nlparams struct is passed in, should be cleaned up since redundant
* `alpha::float64`: optional newmark beta alpha parameter,If TD, use 0.25
* `gamma::float64`: optional newmark beta gamma parameter, if static, use 0. If hydro, use 1.0
* `pBC::Array{<:float}`: Nx3 array consisting of node, local dof, specified displacement value for the boundary condition
* `nodalTerms`: Concentrated nodal terms, should be replaced with the nodal input data array and the calc done internally
* `iterationType::string`: FEA displacement update calculation, Newton Raphson "NR", Direct Iteration "DI"
* `adaptiveLoadSteppingFlag`: Unused, should be removed
* `tolerance::float`: FEA total mesh unsteady analysis convergence tolerance for a timestep
* `maxIterations`: FEA total mesh unsteady analysis convergence max iterations for a timestep
* `maxNumLoadSteps`: used in static (steady state) analysis
* `minLoadStepDelta`: used in static (steady state) analysis
* `minLoadStep`: used in static (steady state) analysis
* `prescribedLoadStep`: used in static (steady state) analysis
* `predef::Bool`: will update the elStorage array if passed into Unsteady() with the nonlinear strain stiffening, to be used for subsequent analyses
* `return_all_reaction_forces::Bool` = return a 2D array of all reaction forces at all timesteps as opposed to the reaction force at a single point.

# Outputs:
* `none`:

"""
function FEAModel(;analysisType = "TNB",
    initCond = [],
    aeroElasticOn = false,
    guessFreq = 0.0,
    airDensity=1.2041,
    gravityOn = true,
    nlOn = false,
    spinUpOn = false,
    outFilename = "none",
    RayleighAlpha = 0.0,
    RayleighBeta = 0.0,
    elementOrder = 1,
    joint = [0,0],
    platformTurbineConnectionNodeNumber = 1,
    jointTransform = 0.0,
    reducedDOFList = zeros(Int,2),
    numDOFPerNode = 6,
    numNodes = 0,
    numModes = 20,
    nlParams = 0,
    alpha = 0.5,
    gamma = 0.5,
    pBC = 0,
    nodalTerms = 0.0,
    iterationType = "NR",
    adaptiveLoadSteppingFlag = true,
    tolerance = 1.0000e-06,
    maxIterations = 50,
    maxNumLoadSteps = 20,
    minLoadStepDelta = 0.0500,
    minLoadStep = 0.0500,
    prescribedLoadStep = 0.0,
    predef = false,
    return_all_reaction_forces = false)

    if jointTransform==0.0
        jointTransform, reducedDOFList = GyricFEA.createJointTransform(joint,numNodes,numDOFPerNode) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    end
    if nodalTerms == 0.0
        # nodalTerms = GyricFEA.readNodalTerms() #Fill in the data structure with nothing
        nodalTerms = applyConcentratedTerms(numNodes, numDOFPerNode)
    end

    if pBC!=0
        BC = makeBCdata(pBC,numNodes,numDOFPerNode,reducedDOFList,jointTransform)
    else
        BC = GyricFEA.BC_struct(0,
        0,
        0,
        0,
        0,
        0,
        0)
    end

    if nlParams==0
        nlParams = NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
        maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep,predef,
        alpha,gamma)
    end

    return FEAModel(analysisType,initCond,aeroElasticOn,guessFreq,airDensity,
    gravityOn,nlOn,spinUpOn,outFilename,RayleighAlpha,RayleighBeta,elementOrder,joint,
    platformTurbineConnectionNodeNumber,jointTransform,reducedDOFList,nlParams,BC,nodalTerms,numModes,return_all_reaction_forces)
end

"""
    NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

See ?FEAModel
"""
mutable struct NlParams
    iterationType
    adaptiveLoadSteppingFlag
    tolerance
    maxIterations
    maxNumLoadSteps
    minLoadStepDelta
    minLoadStep
    prescribedLoadStep
    predef
    alpha
    gamma
end
#Convenience function
NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,maxIterations,maxNumLoadSteps,
minLoadStepDelta,minLoadStep,prescribedLoadStep) = NlParams(iterationType,
adaptiveLoadSteppingFlag,tolerance,maxIterations,maxNumLoadSteps,minLoadStepDelta,
minLoadStep,prescribedLoadStep,false)

"""
Internal, Timoshenko element matrices
"""
mutable struct ElStorage
      K11
      K12
      K13
      K14
      K15
      K16
      K22
      K23
      K24
      K25
      K26
      K33
      K34
      K35
      K36
      K44
      K45
      K46
      K55
      K56
      K66
      M11
      M15
      M16
      M22
      M24
      M33
      M34
      M44
      M55
      M56
      M66
      S11
      S12
      S13
      S15
      S16
      S22
      S23
      S25
      S26
      S33
      S35
      S36
      S55
      S56
      S66
      S14_1
      S14_2
      S24_1
      S24_2
      S34_1
      S34_2
      S45_1
      S45_2
      S46_1
      S46_2
      S44_1
      S44_2
      S44_3
      C12
      C13
      C23
      C24
      C25
      C26
      C34
      C35
      C36
      C14_1
      C14_2
      C45_1
      C45_2
      C46_1
      C46_2
      mel
      moiel
      xmel
      K21NLpredef
      K12NLpredef
      K31NLpredef
      K13NLpredef
      K22NLpredef
      K23NLpredef
      K33NLpredef
end

"""
Internal, time integration terms
"""
mutable struct TimeInt
    delta_t
    a1
    a2
    a3
    a4
    a5
    a6
    a7
    a8
end

"""
Internal, inputs to Timoshenko element, inputs given by FEAModel struct and mesh
"""
mutable struct ElInput
    elementOrder
    modalFlag #not used in timoshenko input?
    timeInt
    xloc
    sectionProps
    sweepAngle
    coneAngle
    rollAngle
    aeroSweepAngle
    iterationType
    useDisp
    preStress
    aeroElasticOn
    aeroForceOn
    loadStepPrev #not used in timoshenko input?
    loadStep
    maxNumLoadSteps #not used in timoshenko input?
    MAXIT #not used in timoshenko input?
    tolerance #not used in timoshenko input?
    analysisType
    disp
    dispdot
    dispddot
    displ_iter
    concMass
    concStiff
    concDamp
    concLoad
    dispm1
    x
    y
    z
    gravityOn
    RayleighAlpha
    RayleighBeta
    accelVec
    omegaVec
    omegaDotVec
    Omega
    OmegaDot
    CN2H
    airDensity
    freq
    firstIteration
end

"""
Internal, timoshenko element output matrices
"""
mutable struct ElOutput
    FhatLessConc
    Ke
    Fe
    Me
    Ce
end

"""
    ElStrain(eps_xx_0,eps_xx_z,eps_xx_y,gam_xz_0,gam_xz_y,gam_xy_0,gam_xy_z)

Struct containing element straing

# Inputs
* `epsilon_x::float`: epsilon_x strain in the x direction
* `epsilon_y::float`: epsilon_y strain in the y direction
* `epsilon_z::float`: epsilon_z strain in the z direction
* `kappa_x::float`: kappa_x curvature in the x direction
* `kappa_y::float`: kappa_y curvature in the y direction
* `kappa_z::float`: kappa_z curvature in the z direction

# Outputs:
* `none`:

"""
mutable struct ElStrain
    epsilon_x
    epsilon_y
    epsilon_z
    kappa_x
    kappa_y
    kappa_z
end

# TODO: get rid of these 3 or four line structs and replace with data arrays
"""
    DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1)

Internal, displacement, velocity, and acceleration for each element

# Inputs
* `elStrain`: Not used, should be removed from this struct
* `displ_sp1::Array{<:float}: displacement position for each dof
* `displddot_sp1::Array{<:float}: displacement acceleration for each dof
* `displdot_sp1::Array{<:float}: displacement velocity for each dof
"""
mutable struct DispOut
    elStrain
    displ_sp1
    displddot_sp1
    displdot_sp1
    eta_sp1
    etadot_sp1
    etaddot_sp1
end
DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1) = DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1,0.0,0.0,0.0)

"""
Internal, displacement, velocity, and acceleration for each element
"""
mutable struct DispData
    displ_s
    displdot_s
    displddot_s
    displ_sm1
    eta_s
    etadot_s
    etaddot_s
end
DispData(displ_s,displdot_s,displddot_s,displ_sm1) = DispData(displ_s,displdot_s,displddot_s,displ_sm1,0.0,0.0,0.0)

"""
    Mesh(nodeNum,numEl,numNodes,x,y,z,elNum,conn,type,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

Struct with mesh definition

# Inputs
* `nodeNum::Array{<:int}`: Number mapping of nodes (typically 1:Nnodes)
* `numEl::int`: total number of elements
* `numNodes::int`: total number of nodes
* `x::Array{<:float}`: Nodal x position
* `y::Array{<:float}`: Nodal y position
* `z::Array{<:float}`: Nodal z position
* `elNum::Array{<:int}`: Number mapping of elements (typically 1:Nelements)
* `conn::Array{<:int}`: Nelemx2 connectivity between nodes, gaps between joints (which are defined in the joints)
* `type::Array{<:int}`: 0-blade 1-tower 2-strut
* `meshSeg::Array{<:int}`: number of nodes within each segment, with segments consisting of tower, blade 1 2 etc, struts
* `structuralSpanLocNorm::Array{<:float}`: Should be named heigh loc norm - unitized position along the blade height, used for aeroload mapping
* `structuralNodeNumbers::Array{<:int}`: Node numbers associated with blades for aero loads mapping
* `structuralElNumbers::Array{<:int}`: Element numbers associated with blades for aero loads mapping

# Outputs:
* `none`:

"""
mutable struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
    type
    meshSeg
    structuralSpanLocNorm
    structuralNodeNumbers
    structuralElNumbers
end

"""
    Ort(Psi_d,Theta_d,Twist_d,Length,elNum,Offset)

Struct with element orientation

# Inputs
* `Psi_d::Array{<:float}`: length NumEl, element rotation about 3 in global FOR (deg) These angles are used to transform from the global coordinate frame to the local element/joint frame via a 3-2 Euler rotation sequence.
* `Theta_d::Array{<:float}`: length NumEl, element rotation about 2 (deg)
* `Twist_d::Array{<:float}`: length NumEl, element twist (deg)
* `Length::Array{<:float}`: length NumEl, element length (m)
* `elNum::Array{<:float}`: Element number the other arrays are associated with
* `Offset::Array{<:float}`: hub frame coordinate of node 1 of the element

# Outputs:
* `none`:

"""
mutable struct Ort
    Psi_d
    Theta_d
    Twist_d
    Length
    elNum
    Offset
end

"""
Internal, boundary condition data, see ?FEAModel for pBC
"""
mutable struct BC_struct
    numpBC
    pBC
    numsBC
    nummBC
    isConstrained
    map
    redVectorMap
end

"""
    SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset,xaf,yaf)

Struct with element sectional properties, each component is a 1x2 array with distributed properties

# Inputs
* `ac::Array{<:float}`: aerodynamic center, used in flutter approximation
* `twist::Array{<:float}`: element twist (rad)
* `rhoA::Array{<:float}`: rho * A in standard SI units
* `EIyy::Array{<:float}`: E * Iyy
* `EIzz::Array{<:float}`: E * Izz
* `GJ::Array{<:float}`: G * J
* `EA::Array{<:float}`: E * A
* `rhoIyy::Array{<:float}`: rho * Iyy
* `rhoIzz::Array{<:float}`: rho * Izz
* `rhoJ::Array{<:float}`: rho * J
* `zcm::Array{<:float}`: z location of center of mass
* `ycm::Array{<:float}`: y location of center of mass
* `a::Array{<:float}`: possibly lift slope
* `EIyz::Array{<:float}`: E * Iyz
* `alpha1::Array{<:float}`: #This is always 0 in the element file, and it is unclear what it is used for since I can't find it being used in the code
* `alpha2::Array{<:float}`: doesn't appear to be used
* `alpha3::Array{<:float}`: doesn't appear to be used
* `alpha4::Array{<:float}`: doesn't appear to be used
* `alpha5::Array{<:float}`: doesn't appear to be used
* `alpha6::Array{<:float}`: doesn't appear to be used
* `rhoIyz::Array{<:float}`: rho * Iyz
* `b::Array{<:float}`: used in flutter approximation, possibly a chord or thickness value
* `a0::Array{<:float}`: zero lift angle of attack, used in flutter approximation
* `aeroCenterOffset::Array{<:float}`: doesn't appear to be used
* `xaf::Array{<:float}`: x airfoil coordinates (to scale)
* `yaf::Array{<:float}`: y airfoil coordinates (to scale)

# Outputs:
* `none`:

"""
mutable struct SectionPropsArray
    ac
    twist
    rhoA
    EIyy
    EIzz
    GJ
    EA
    rhoIyy
    rhoIzz
    rhoJ
    zcm
    ycm
    a
    EIyz
    alpha1
    alpha2
    alpha3
    alpha4
    alpha5
    alpha6
    rhoIyz
    b
    a0
    aeroCenterOffset
    xaf
    yaf
end
SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset) = SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset,nothing,nothing) #convenience function

"""
Internal, see ?Ort and ?SectionPropsArray
"""
mutable struct El
    props
    elLen
    psi
    theta
    roll
    rotationalEffects
end

"""
Internal, see ?FEAModel for NodalTerms
"""
mutable struct NodalTerms
    concLoad
    concStiff
    concMass
    concDamp
end

"""
Internal, NodalTerms node number, local dof (diagonal), and value
"""
mutable struct ConcNDL1D
    nodeNum
    dof
    val
end

"""
Internal, NodalTerms node number, local dof1, local dof2, and value
"""
mutable struct ConcNDL2D
    nodeNum
    dof1
    dof2
    val
end

"""
Internal, ROM data
"""
mutable struct ROM
    Kr
    Mr
    Cr
    CrModal
    Fr
    Phi
    invPhi
    SrOx2
    SrOy2
    SrOz2
    SrOxOy
    SrOyOz
    SrOxOz
    FrOx2
    FrOy2
    FrOz2
    FrOxOy
    FrOyOz
    FrOxOz
    FrOxdot
    FrOydot
    FrOzdot
    FrAx
    FrAy
    FrAz
    GrOx
    GrOy
    GrOz
    HrOx
    HrOy
    HrOz
end

ROM(Kr,Cr,Fr) = ROM(Kr,0.0,Cr,0.0,Fr,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ROM(Kr,Mr,Cr,CrModal,Fr,Phi,invPhi) = ROM(Kr,Mr,Cr,CrModal,Fr,Phi,invPhi,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
