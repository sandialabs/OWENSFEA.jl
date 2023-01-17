"""
    calculateTimoshenkoElementInitialRun(elementOrder,modalFlag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,x,y,z,concMassFlag,concMass,Omega)

Internal, see ?initialElementCalculations, performs initial element calculations and stores them for later use and efficiency gains.
"""
function calculateTimoshenkoElementInitialRun(elementOrder,modalFlag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,x,y,z,concMassFlag,concMass,Omega)

    numGP = 4

    #calculate quad points
    xi,weight = getGP(numGP)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    K11 = zeros(numNodesPerEl,2)
    K12 = zero(K11)
    K13 = zero(K11)
    K14 = zero(K11)
    K15 = zero(K11)
    K16 = zero(K11)
    K22 = zero(K11)
    K23 = zero(K11)
    K24 = zero(K11)
    K25 = zero(K11)
    K26 = zero(K11)
    K33 = zero(K11)
    K34 = zero(K11)
    K35 = zero(K11)
    K36 = zero(K11)
    K44 = zero(K11)
    K45 = zero(K11)
    K46 = zero(K11)
    K55 = zero(K11)
    K56 = zero(K11)
    K66 = zero(K11)

    S11 = zero(K11)
    S12 = zero(K11)
    S13 = zero(K11)
    S14_1 = zero(K11)
    S14_2 = zero(K11)
    S15 = zero(K11)
    S16 = zero(K11)
    S22 = zero(K11)
    S23 = zero(K11)
    S24_1 = zero(K11)
    S24_2 = zero(K11)
    S25 = zero(K11)
    S26 = zero(K11)
    S33 = zero(K11)
    S34_1 = zero(K11)
    S34_2 = zero(K11)
    S35 = zero(K11)
    S36 = zero(K11)
    S44_1 = zero(K11)
    S44_2 = zero(K11)
    S44_3 = zero(K11)
    S45_1 = zero(K11)
    S45_2 = zero(K11)
    S46_1 = zero(K11)
    S46_2 = zero(K11)
    S55 = zero(K11)
    S56 = zero(K11)
    S66 = zero(K11)

    #     F1 = zeros(numNodesPerEl,1)
    #     F3 = F1
    #     F2 = F1
    #     F4 = F1
    #     F5 = F1
    #     F6 = F1


    M11 = zero(K11)
    M15 = zero(K11)
    M16 = zero(K11)
    M22 = zero(K11)
    M24 = zero(K11)
    M33 = zero(K11)
    M34 = zero(K11)
    M44 = zero(K11)
    M55 = zero(K11)
    M56 = zero(K11)
    M66 = zero(K11)

    C12 = zero(K12)
    C13 = zero(K13)
    C14_1 = zero(K14)
    C14_2 = zero(K14)
    C23 = zero(K23)
    C24 = zero(K24)
    C34 = zero(K34)
    C25 = zero(K11)
    C26 = zero(K11)
    C35 = zero(K11)
    C36 = zero(K11)
    C45_1 = zero(K11)
    C45_2 = zero(K11)
    C46_1 = zero(K11)
    C46_2 = zero(K11)

    elementMass = 0.0
    elementItens = zeros(3,3)
    elxm = zeros(3)

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg_d = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg_d.*pi/180.0)

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = N
        p_N1_x = p_N_x
        N2 = N
        p_N2_x = p_N_x
        N3 = N
        p_N3_x = p_N_x
        N4 = N
        p_N4_x = p_N_x
        N5 = N
        p_N5_x = p_N_x
        N6 = N
        p_N6_x = p_N_x
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        EA   = interpolateVal(sectionProps.EA,N) #struct stiffness terms
        EIyy = interpolateVal(sectionProps.EIyy,N)
        EIzz = interpolateVal(sectionProps.EIzz,N)
        GJ   = interpolateVal(sectionProps.GJ,N)
        EIyz = interpolateVal(sectionProps.EIyz,N)

        couple16 = 0.0   #int(Ey dA)    v bend - extension
        couple15 = 0.0   #int(Ez dA)    w bend - extension
        couple45 = 0.0   #int(Gz dA)    w bend - twist
        couple46 = 0.0   #int(Gz dA)    v bend - twist
        couple14 = 0.0   # extension twist
        couple34 = couple46
        couple24 = couple45

        rhoA   = interpolateVal(sectionProps.rhoA,N) #struct mass terms
        rhoIyy = interpolateVal(sectionProps.rhoIyy,N)
        rhoIzz = interpolateVal(sectionProps.rhoIzz,N)
        rhoJ   = interpolateVal(sectionProps.rhoJ,N)
        rhoIyz = interpolateVal(sectionProps.rhoIyz,N)

        vprime = 0.0 #set to zero to deactivate nonlinearites from initial element calculations
        wprime = 0.0

        ycm = interpolateVal(sectionProps.ycm,N)
        zcm = interpolateVal(sectionProps.zcm,N)

        xgp      = interpolateVal(x,N1)
        ygp      = interpolateVal(y,N1)
        zgp      = interpolateVal(z,N1)

        #.... end interpolate value at quad points ........

        #adjust moments of inertia for offsets
        rhoIyy = rhoIyy + rhoA*zcm^2
        rhoIzz = rhoIzz + rhoA*ycm^2
        rhoIyz = rhoIyz + rhoA*ycm*zcm
        rhoJ   = rhoJ + rhoA*(ycm^2 + zcm^2)

        #Calculate strutural stiffness sub matrices
        calculateElement1!(EA,integrationFactor,p_N1_x,p_N1_x,K11)
        calculateElement1!(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12)
        calculateElement1!(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13)
        calculateElement1!(couple14,integrationFactor,p_N1_x,p_N4_x,K14)
        calculateElement1!(couple15,integrationFactor,p_N1_x,p_N5_x,K15)
        calculateElement1!(-couple16,integrationFactor,p_N1_x,p_N6_x,K16)
        calculateElement1!(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22)
        calculateElement1!(-couple24,integrationFactor,p_N2_x,p_N4_x,K24)
        calculateElement1!(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33)
        calculateElement1!(couple34,integrationFactor,p_N3_x,p_N4_x,K34)
        calculateElement1!(GJ,integrationFactor,p_N4_x,p_N4_x,K44)
        calculateElement1!(couple45,integrationFactor,p_N4_x,N5,K45)
        calculateElement1!(couple46,integrationFactor,p_N4_x,N6,K46)
        calculateElement1!(EIyy,integrationFactor,p_N5_x,p_N5_x,K55)
        calculateElement1!(-EIyz,integrationFactor,p_N5_x,p_N6_x,K56)
        calculateElement1!(EIzz,integrationFactor,p_N6_x,p_N6_x,K66)

        #Calculate structural mass sub matrices
        calculateElement1!(rhoA,integrationFactor,N1,N1,M11)
        calculateElement1!(rhoA*zcm,integrationFactor,N1,N5,M15)
        calculateElement1!(-rhoA*ycm,integrationFactor,N1,N6,M16)
        calculateElement1!(rhoA,integrationFactor,N2,N2,M22)
        calculateElement1!(-rhoA*zcm,integrationFactor,N2,N4,M24)
        calculateElement1!(rhoA,integrationFactor,N3,N3,M33)
        calculateElement1!(rhoA*ycm,integrationFactor,N3,N4,M34)
        calculateElement1!(rhoJ,integrationFactor,N4,N4,M44)
        calculateElement1!(rhoIyy,integrationFactor,N5,N5,M55)
        calculateElement1!(-rhoIyz,integrationFactor,N5,N6,M56)
        calculateElement1!(rhoIzz,integrationFactor,N6,N6,M66)

        #Calculate Centrifugal load vector and gravity load vector
        #eventually incorporate lambda into gp level to account for variable
        #twist

        O1 = 1 #these are set to unity to get coefficients for omega components
        O2 = 1
        O3 = 1

        posLocal = lambda[1:3,1:3]*[xgp, ygp, zgp]
        xbarlocal = posLocal[1]
        ybarlocal = posLocal[2]
        zbarlocal = posLocal[3]

        #        g=9.81 #gravitational acceleration [m/s^2]
        #        a_x = 0 #acceleration of body in x and y (hardwired to zero for now)
        #        a_y = 0
        #        a_z = -g
        #        fx = rhoA*a_x #let these loads be defined in the inertial frame
        #        fy = rhoA*a_y
        #        fz = rhoA*a_z
        #        rvec = [ 0 ycm zcm]
        #
        #        fi_hub = CN2H*[fxfyfz]
        #
        #        disLoadgpLocal = lambda(1:3,1:3)*fi_hub
        #        disMomentgp = cross(rvec,disLoadgpLocal)

        #        f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal) - disLoadgpLocal[1]    #omega dot loading not
        #        calculateVec1!(f1,integrationFactor,N1,F1)
        #        f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2) - disLoadgpLocal[2]
        #        calculateVec1!(f2,integrationFactor,N2,F2)
        #        f3 = rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal) - disLoadgpLocal[3]
        #        calculateVec1!(f3,integrationFactor,N3,F3)
        #        f4 = rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2))...
        #                   + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3)) - disMomentgp[1]
        #        calculateVec1!(f4,integrationFactor,N4,F4)
        #        f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3) - disMomentgp[2]
        #        calculateVec1!(f5,integrationFactor,N5,F5)
        #        f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2))) - disMomentgp[3]
        #        calculateVec1!(f6,integrationFactor,N6,F6)
        #
        #Gyric matrix (Coriolis)
        calculateElement1!(-2*rhoA*O3,integrationFactor,N1,N2,C12)
        calculateElement1!(2*rhoA*O2,integrationFactor,N1,N3,C13)

        calculateElement1!(2*rhoA*(ycm*O2),integrationFactor,N1,N4,C14_1)
        calculateElement1!(2*rhoA*(zcm*O3),integrationFactor,N1,N4,C14_2)

        calculateElement1!(-2*rhoA*O1,integrationFactor,N2,N3,C23)
        calculateElement1!(-2*rhoA*ycm*O1,integrationFactor,N2,N4,C24)
        calculateElement1!(2*rhoA*zcm*O3,integrationFactor,N2,N5,C25)
        calculateElement1!(-2*rhoA*ycm*O3,integrationFactor,N2,N6,C26)
        calculateElement1!(-2*rhoA*zcm*O1,integrationFactor,N3,N4,C34)
        calculateElement1!(-2*rhoA*zcm*O2,integrationFactor,N3,N5,C35)
        calculateElement1!(2*rhoA*ycm*O2,integrationFactor,N3,N6,C36)

        calculateElement1!(-2*(rhoIyy*O3),integrationFactor,N4,N5,C45_1)
        calculateElement1!(-2*(rhoIyz*O2),integrationFactor,N4,N5,C45_2)

        calculateElement1!(2*(rhoIzz*O2),integrationFactor,N4,N6,C46_1)
        calculateElement1!(2*(rhoIyz*O3),integrationFactor,N4,N6,C46_2)

        #Spin softening matrix
        calculateElement1!(-rhoA*(O2^2+O3^2),integrationFactor,N1,N1,S11)
        calculateElement1!(rhoA*O1*O2,integrationFactor,N1,N2,S12)
        calculateElement1!(rhoA*O1*O3,integrationFactor,N1,N3,S13)

        calculateElement1!(rhoA*(ycm*O1*O3),integrationFactor,N1,N4,S14_1)
        calculateElement1!(rhoA*(-zcm*O1*O2),integrationFactor,N1,N4,S14_2)

        calculateElement1!(-rhoA*zcm*(O2^2+O3^2),integrationFactor,N1,N5,S15)
        calculateElement1!(rhoA*ycm*(O2^2+O3^2),integrationFactor,N1,N6,S16)
        calculateElement1!(-rhoA*(O1^2+O3^2),integrationFactor,N2,N2,S22)
        calculateElement1!(rhoA*O2*O3,integrationFactor,N2,N3,S23)

        calculateElement1!(rhoA*zcm*(O1^2+O3^2),integrationFactor,N2,N4,S24_1)
        calculateElement1!(rhoA*ycm*O2*O3,integrationFactor,N2,N4,S24_2)

        calculateElement1!(rhoA*zcm*O1*O2,integrationFactor,N2,N5,S25)
        calculateElement1!(-rhoA*ycm*O1*O2,integrationFactor,N2,N6,S26)
        calculateElement1!(-rhoA*(O1^2+O2^2),integrationFactor,N3,N3,S33)

        calculateElement1!(-rhoA*(ycm*(O1^2+O2^2)),integrationFactor,N3,N4,S34_1)
        calculateElement1!(-rhoA*(zcm*O2*O3),integrationFactor,N3,N4,S34_2)


        calculateElement1!(rhoA*zcm*O1*O3,integrationFactor,N3,N5,S35)
        calculateElement1!(-rhoA*ycm*O1*O3,integrationFactor,N3,N6,S36)

        calculateElement1!(-(rhoIyy*(O1^2+O3^2)),integrationFactor,N4,N4,S44_1)
        calculateElement1!(-(rhoIzz*(O1^2+O2^2)),integrationFactor,N4,N4,S44_2)
        calculateElement1!(-(2*rhoIyz*O2*O3),integrationFactor,N4,N4,S44_3)

        calculateElement1!(rhoIyz*O1*O3,integrationFactor,N4,N5,S45_1)
        calculateElement1!(-rhoIyy*O1*O2,integrationFactor,N4,N5,S45_2)

        calculateElement1!(rhoIyz*O1*O2,integrationFactor,N4,N6,S46_1)
        calculateElement1!(-rhoIzz*O1*O3,integrationFactor,N4,N6,S46_2)


        calculateElement1!(-rhoIyy*(O2^2+O3^2),integrationFactor,N5,N5,S55)
        calculateElement1!(rhoIyz*(O2^2+O3^2),integrationFactor,N5,N6,S56)
        calculateElement1!(-rhoIzz*(O2^2+O3^2),integrationFactor,N6,N6,S66)

        elementMass,elementItens,elxm = calculateElementMass(rhoA,rhoIyy,rhoIzz,rhoIyz,rhoJ,ycm,zcm,xbarlocal,ybarlocal,zbarlocal,integrationFactor,elementMass,elementItens,elxm)

    end

    #######################################
    #Reduced integration loop
    numGP = 1
    xi,weight = getGP(numGP)

    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N5 = N
        N6 = N
        p_N2_x = p_N_x
        p_N3_x = p_N_x
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        EA   = interpolateVal(sectionProps.EA,N) #struct stiffness terms
        GA = EA/2.6*5/6
        #.... end interpolate value at quad points ........

        #Calculate strutural stiffness sub matrices
        calculateElement1!(GA,integrationFactor,p_N2_x,p_N2_x,K22)
        calculateElement1!(-GA,integrationFactor,p_N2_x,N6,K26)
        calculateElement1!(GA,integrationFactor,p_N3_x,p_N3_x,K33)
        calculateElement1!(GA,integrationFactor,p_N3_x,N5,K35)
        calculateElement1!(GA,integrationFactor,N5,N5,K55)
        calculateElement1!(GA,integrationFactor,N6,N6,K66)

    end

    lamSlim = lambda[1:3,1:3]
    lamSlimTran = lamSlim'
    elementMOI = lamSlimTran*elementItens*lamSlim[1:3,1:3]
    elxm = lamSlimTran*elxm
    #
    #
    if (concMassFlag)
        #modify element mass, moi, and xm to account for concentrated terms
        elementMass = elementMass + concMass[1,1] + concMass[1,7] #TODO: why not the other terms? + sum(concMass[1,:])

        elementMOI[1,1] = elementMOI[1,1] + concMass[1,1]*(y[1]^2 + z[1]^2)+ concMass[1,7]*(y[2]^2 + z[2]^2) + concMass[4,4] + concMass[4,10]
        elementMOI[2,2] = elementMOI[2,2] + concMass[1,1]*(x[1]^2 + z[1]^2)+ concMass[1,7]*(x[2]^2 + z[2]^2) + concMass[5,5] + concMass[5,11]
        elementMOI[3,3] = elementMOI[3,3] + concMass[1,1]*(x[1]^2 + y[1]^2)+ concMass[1,7]*(x[2]^2 + y[2]^2) + concMass[6,6] + concMass[6,12]
        elementMOI[1,2] = elementMOI[1,2] - concMass[1,1]*x[1]*y[1] - concMass[1,7]*x[2]*y[2]
        elementMOI[1,3] = elementMOI[1,3] - concMass[1,1]*x[1]*z[1] - concMass[1,7]*x[2]*z[2]
        elementMOI[2,1] = elementMOI[2,1] - concMass[1,1]*x[1]*y[1] - concMass[1,7]*x[2]*y[2]
        elementMOI[2,3] = elementMOI[2,3] - concMass[1,1]*y[1]*z[1] - concMass[1,7]*y[2]*z[2]
        elementMOI[3,1] = elementMOI[3,1] - concMass[1,1]*x[1]*z[1] - concMass[1,7]*x[2]*z[2]
        elementMOI[3,2] = elementMOI[3,2] - concMass[1,1]*y[1]*z[1] - concMass[1,7]*y[2]*z[2]

        elxm[1] = elxm[1] + concMass[1,1]*x[1] + concMass[1,7]*x[2]
        elxm[2] = elxm[2] + concMass[1,1]*y[1] + concMass[1,7]*y[2]
        elxm[3] = elxm[3] + concMass[1,1]*z[1] + concMass[1,7]*z[2]
    end

    #store element properties
    return ElStorage(K11,     #Store structural stiffness "K" into elementStorage
    K12,K13,K14,K15,K16,K22,K23,K24,K25,K26,K33,K34,K35,K36,K44,K45,K46,K55,K56,K66,
    M11, #Store structural stiffness "M" into elementStorage
    M15,M16,M22,M24,M33,M34,M44,M55,M56,M66,
    0.5*S11, #Store spin softening coefficient "S" into element storage
    S12,S13,0.5*S15,0.5*S16,0.5*S22,S23,S25,S26,0.5*S33,S35,S36,0.5*S55,0.5*S56,0.5*S66,S14_1,S14_2,S24_1,S24_2,S34_1,S34_2,S45_1,S45_2,S46_1,S46_2,S44_1,S44_2,S44_3,
    C12, #Store coriolis coefficient "C" into element sotrage
    C13,C23,C24,C25,C26,C34,C35,C36,C14_1,C14_2,C45_1,C45_2,C46_1,C46_2,
    elementMass,
    elementMOI,
    elxm,zeros(2,2),zeros(2,2),zeros(2,2),zeros(2,2),zeros(2,2),zeros(2,2),zeros(2,2))
end


"""
    calculateTimoshenkoElementNL(input,elStorage;predef=nothing)

Internal, performs nonlinear element calculations.

#Inputs
* `input::ElInput`: see ?ElInput
* `elStorage::ElStorage`: see ?ElStorage
* `predef::Bool`: optional, if true, mutates ElStorage to include the nonlinear strain stiffening

#Outputs
* `ElOutput`: see ?ElOutput

"""
function calculateTimoshenkoElementNL(input,elStorage;predef=nothing)

    ###-------- assign input block ----------------
    elementOrder   = input.elementOrder
    x              = input.x
    y              = input.y
    z              = input.z
    xloc           = input.xloc
    disp           = input.disp
    sectionProps   = input.sectionProps
    sweepAngle     = input.sweepAngle
    coneAngle      = input.coneAngle
    rollAngle      = input.rollAngle
    Omega          = input.Omega
    OmegaDot       = input.OmegaDot
    concStiff      = input.concStiff
    concMass       = input.concMass
    concDamp       = input.concDamp
    concLoad       = input.concLoad
    useDisp        = input.useDisp
    preStress      = input.preStress
    aeroElasticOn  = input.aeroElasticOn
    aeroForceOn    = input.aeroForceOn
    iterationType  = input.iterationType
    disp_iter      = input.displ_iter
    omegaVec       = input.omegaVec
    CN2H           = input.CN2H
    omegaDotVec    = input.omegaDotVec
    timeInt        = input.timeInt
    analysisType = input.analysisType
    dispm1 = zeros(12) #declare type
    dispdot = zeros(12)#declare type
    dispddot = zeros(12)#declare type
    accelVec       = input.accelVec

    #options for Dean integrator
    if analysisType=="TD"
        dispm1         = input.dispm1
    elseif analysisType=="TNB"#options for newmark beta integrator
        dispdot      = input.dispdot
        dispddot     = input.dispddot
    end

    #--------------------------------------------
    #setting for modal analysis flag
    if analysisType=="M"
        disp_iter=copy(disp)
    end

    #setting for initial reduced order feamodel calculations
    if analysisType=="ROM"
        disp_iter=copy(disp)
        omegaVec = input.omegaVec
        omegaDotVec = input.omegaDotVec
    end

    #settings if aeroelastic analysis is active
    if aeroElasticOn
        freq = input.freq
    else
        freq = 0.0 #Not used, but must be declared
    end


    numGP = 4   #number of gauss points for full integration
    numGPRI = 1 #number of gauss points for reduced integration
    #calculate quad points
    xi,weight = getGP(numGP)
    xiRI,weightRI = getGP(numGPRI)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    F1 = zeros(numNodesPerEl)
    F3 = zero(F1)
    F2 = zero(F1)
    F4 = zero(F1)
    F5 = zero(F1)
    F6 = zero(F1)

    SS22 = zeros(numNodesPerEl,numNodesPerEl) #initialize pre-stress (stress stiffening matrices)
    SS33 = zero(SS22)

    #initialize nonlinear element matrices, only used if (useDisp)
    K12NL = zeros(numNodesPerEl,numNodesPerEl)
    K13NL = zero(K12NL)
    K22NL = zero(K12NL)
    K23NL = zero(K12NL)
    K33NL = zero(K12NL)
    K22NLhat = zero(K22NL)
    K33NLhat = zero(K33NL)
    K23NLhat = zero(K23NL)


    #initialize aeroelastic matrices only used if aeroElasticOn, but must declare type
    K33Aero = zeros(numNodesPerEl,numNodesPerEl)
    K34Aero = zero(K33Aero)
    C33Aero = zero(K33Aero)
    C34Aero = zero(K33Aero)
    K43Aero = zero(K33Aero)
    K44Aero = zero(K33Aero)
    C43Aero = zero(K33Aero)
    C44Aero = zero(K33Aero)

    C33 = zeros(numNodesPerEl,numNodesPerEl)
    C44 = zero(C33)

    #Convert frequencies from Hz to radians
    Omega = 2*pi*Omega
    OmegaDot = 2*pi*OmegaDot

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node

    twistAvg_d = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg_d.*pi/180.0)
    lambdaSlim = lambda[1:3,1:3]

    dispLocal = lambda*disp_iter

    uNode = [dispLocal[1] dispLocal[7]]
    vNode = [dispLocal[2] dispLocal[8]]
    wNode = [dispLocal[3] dispLocal[9]]
    #     theta_xNode = [dispLocal(4)  dispLocal(10)]
    #     theta_yNode = [dispLocal(5)  dispLocal(11)]
    #     theta_zNode = [dispLocal(6)  dispLocal(12)]


    omega_x=omegaVec[1]
    omega_y=omegaVec[2]
    omega_z = omegaVec[3] + Omega
    omegaDot_x=omegaDotVec[1]
    omegaDot_y=omegaDotVec[2]
    omegaDot_z = omegaDotVec[3] + OmegaDot
    Ohub = [omega_x;omega_y;omega_z]
    ODotHub = [omegaDot_x;omegaDot_y;omegaDot_z]
    Oel = lambdaSlim*Ohub
    ODotel = lambdaSlim*ODotHub
    O1 = Oel[1]
    O2 = Oel[2]
    O3 = Oel[3]

    O1dot = ODotel[1]
    O2dot = ODotel[2]
    O3dot = ODotel[3]

    if eltype(input.gravityOn) == Bool && input.gravityOn == true
        a_x_n = 0.0 #accelerations in inertial frame
        a_y_n = 0.0
        a_z_n = 9.81 # gravity
    elseif eltype(input.gravityOn) == Bool && input.gravityOn == false
        a_x_n = 0.0 #accelerations in inertial frame
        a_y_n = 0.0
        a_z_n = 0.0
    end

    if eltype(input.gravityOn) == Float64
        a_x_n = input.gravityOn[1] #accelerations in inertial frame
        a_y_n = input.gravityOn[2]
        a_z_n = input.gravityOn[3]
    end

    a_x = accelVec[1] #acceleration of body in hub frame (from platform rigid body motion)
    a_y = accelVec[2]
    a_z = accelVec[3]

    a_temp = CN2H*[a_x_n; a_y_n; a_z_n]

    a_x = a_x + a_temp[1]
    a_y = a_y + a_temp[2]
    a_z = a_z + a_temp[3]


    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = copy(N);  p_N1_x = copy(p_N_x)
        N2 = copy(N);  p_N2_x = copy(p_N_x)
        N3 = copy(N);  p_N3_x = copy(p_N_x)
        N4 = copy(N);
        N5 = copy(N);
        N6 = copy(N);
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        rhoA = interpolateVal(sectionProps.rhoA,N) #struct mass terms

        # Only used if (useDisp || preStress)
        EA = interpolateVal(sectionProps.EA,N)
        uprime = interpolateVal(uNode,p_N1_x)
        vprime = interpolateVal(vNode,p_N2_x)
        wprime = interpolateVal(wNode,p_N3_x)
        Faxial = EA*(uprime+0.5*vprime^2+0.5*wprime^2)


        #mass center offsets
        ycm = interpolateVal(sectionProps.ycm,N)
        zcm = interpolateVal(sectionProps.zcm,N)

        xgp = interpolateVal(x,N1)
        ygp = interpolateVal(y,N1)
        zgp = interpolateVal(z,N1)

        if aeroElasticOn || aeroForceOn
            #aerodynamic props/data
            airDensity = input.airDensity
            acgp     = interpolateVal(sectionProps.ac,N1)
            a0gp     = interpolateVal([sectionProps.a0],N1)
            bgp      = interpolateVal(sectionProps.b,N1)
            agp      = interpolateVal(sectionProps.a,N1)
            radiusgp = sqrt(xgp^2+ygp^2) # radiusgp = 0 for tower
            Uinfgp   = Omega*radiusgp
            twistgp = interpolateVal(sectionProps.twist,N1)
            #.... end interpolate value at quad points ........
        else
            airDensity = 0.0 #Not used but must declare type
            acgp     = 0.0 #Not used but must declare type
            a0gp     = 0.0 #Not used but must declare type
            bgp      = 0.0 #Not used but must declare type
            agp      = 0.0 #Not used but must declare type
            Uinfgp   = 0.0 #Not used but must declare type
            twistgp = 0.0 #Not used but must declare type
        end

        if aeroElasticOn
            kgp      = freq*bgp/Uinfgp
            Theogp   = calculateTheo(kgp)
        else
            kgp = 0.0 #Not used but must declare type
            Theogp = 0.0 #Not used but must declare type
        end


        #Calculate Centrifugal load vector and gravity load vector
        #eventually incorporate lambda into gp level to account for variable
        #twist
        posLocal = lambdaSlim*[xgp;ygp;zgp]
        xbarlocal = posLocal[1]
        ybarlocal = posLocal[2]
        zbarlocal = posLocal[3]

        fx = rhoA*a_x #let these loads be defined in the inertial frame
        fy = rhoA*a_y
        fz = rhoA*a_z
        rvec = [ 0; ycm; zcm]

        fi_hub = [fx;fy;fz]

        disLoadgpLocal = lambdaSlim*fi_hub
        cpskew= [0 -rvec[3] rvec[2]
                rvec[3] 0 -rvec[1]
                -rvec[2] rvec[1] 0]
        disMomentgp = cpskew*disLoadgpLocal

        if preStress #stress-stiffening/pre-stress calculations
            calculateElement1!(Faxial,integrationFactor,p_N2_x,p_N2_x,SS22)
            calculateElement1!(Faxial,integrationFactor,p_N3_x,p_N3_x,SS33)
        end


        #calculate static aerodynamic load
        sectionAeroLift = 0.0
        sectionAeroMoment = 0.0
        if (aeroForceOn)
            cl = a0gp*twistgp*pi/180.0
            qinf = 0.5*airDensity*Uinfgp^2*(2.0*bgp)
            sectionAeroLift = qinf*cl
            sectionAeroMoment = sectionAeroLift*(acgp+agp)
        end

        #distributed/body force load calculations
        f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal + O3dot*ybarlocal - O2dot*zbarlocal) - disLoadgpLocal[1]
        calculateVec1!(f1,integrationFactor,N1,F1)
        f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2 + O1dot*zbarlocal - O3dot*xbarlocal) - disLoadgpLocal[2]
        calculateVec1!(f2,integrationFactor,N2,F2)
        f3 = sectionAeroLift + rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal + O2dot*xbarlocal - O1dot*ybarlocal) - disLoadgpLocal[3]
        calculateVec1!(f3,integrationFactor,N3,F3)
        f4 = sectionAeroMoment + rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2)) + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3) + ycm*(O2dot*xbarlocal - O1dot*ybarlocal) - zcm*(O1dot*zbarlocal - O3dot*xbarlocal)) - disMomentgp[1]
        calculateVec1!(f4,integrationFactor,N4,F4)
        f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3 - O2dot*zbarlocal + O3dot*ybarlocal) - disMomentgp[2]
        calculateVec1!(f5,integrationFactor,N5,F5)
        f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2)) - O3dot*ybarlocal + O2dot*zbarlocal) - disMomentgp[3]
        calculateVec1!(f6,integrationFactor,N6,F6)

        if aeroElasticOn && (bgp != 0) #aeroelastic calculations
            #This is a real valued aeroelastic representation from
            #Wright and Cooper
            #This version assumes aerodynamic center at quarter chord of
            #airfoil
            Fgp = real(Theogp)
            Ggp = imag(Theogp)
            lcsrat = a0gp/(2*pi)
            Fgp = Fgp*lcsrat
            Ggp = Ggp*lcsrat
            agp = agp/bgp

            if (Uinfgp==0)
                kgp = 1
            end
            lz = -2*pi*(-0.5*kgp^2-Ggp*kgp) #leading negative for difference in unsteady z and w-flap direction
            lzdot = -2*pi*Fgp #leading negative for difference in unsteady z and w-flap direction
            ltheta = -2*pi*(0.5*kgp^2*agp + Fgp - Ggp*kgp*(0.5-agp)) #leading negative for difference in pitch and torsion angles
            if (kgp==0)
                lthetadot = -2*pi*(.5 + Fgp*(.5-agp))
            else
                lthetadot = -2*pi*(.5 + Fgp*(.5-agp) + Ggp/kgp)
            end

            mz = -2*pi*(-0.5*kgp^2*agp-kgp*(agp+acgp)*Ggp) #same as above for lz,lzdot,leheta,lthetadot
            mzdot = -2*pi*(agp+acgp)*Fgp
            mtheta = -2*pi*(0.5*kgp^2*(1.0/8.0+agp^2)+Fgp*(agp+acgp)-kgp*Ggp*(agp+0.5)*(0.5-agp))
            if (kgp==0)
                mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp))
            else
                mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp)+Ggp/kgp*(agp+acgp))
            end

            k33fac = airDensity*Uinfgp^2*lz
            k34fac = airDensity*Uinfgp^2*bgp*ltheta
            k43fac = -airDensity*Uinfgp^2*bgp*mz #leading negative
            k44fac = -airDensity*Uinfgp^2*bgp^2*mtheta #leading negative

            c33fac = airDensity*Uinfgp*bgp*lzdot
            c34fac = airDensity*Uinfgp*bgp^2*lthetadot
            c43fac = -airDensity*Uinfgp*bgp^2*mzdot #leading negative
            c44fac = -airDensity*Uinfgp*bgp^3*mthetadot #leading negative

            calculateElement1!(-k33fac,integrationFactor,N3,N3,K33Aero)
            calculateElement1!(-k34fac,integrationFactor,N3,N4,K34Aero)
            calculateElement1!(-c34fac,integrationFactor,N3,N4,C34Aero)
            calculateElement1!(-c33fac,integrationFactor,N3,N3,C33Aero)

            calculateElement1!(-k43fac,integrationFactor,N4,N3,K43Aero)
            calculateElement1!(-k44fac,integrationFactor,N4,N4,K44Aero)
            calculateElement1!(-c43fac,integrationFactor,N4,N3,C43Aero)
            calculateElement1!(-c44fac,integrationFactor,N4,N4,C44Aero)

        end

    end #END OF INTEGRATION LOOP

    #Integration loop
    for i=1:numGPRI
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xiRI[i],xloc)
        p_N1_x = copy(p_N_x)
        p_N2_x = copy(p_N_x)
        p_N3_x = copy(p_N_x)

        integrationFactor = Jac * weightRI[i]

        #..... interpolate for value at quad point .....

        if useDisp || preStress
            EA   = interpolateVal(sectionProps.EA,N)
            uprime = interpolateVal(uNode,p_N1_x)
            vprime = interpolateVal(vNode,p_N2_x)
            wprime = interpolateVal(wNode,p_N3_x)
        end


        if useDisp
            #nonlinear element calculations
            calculateElement1!(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12NL)
            calculateElement1!(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13NL)
            calculateElement1!(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22NL)
            calculateElement1!(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NL)
            calculateElement1!(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33NL)

            #K12NLhat = K12
            #K13NLhat = K13
            #nonlinear element tangent matrix component calculations
            # T_ij = K_ij + Khat_ij
            if iterationType=="NR"
                calculateElement1!(EA*(uprime + vprime^2 + 0.5*wprime^2),integrationFactor,p_N2_x,p_N2_x,K22NLhat)
                calculateElement1!(EA*(uprime + wprime^2 + 0.5*vprime^2),integrationFactor,p_N3_x,p_N3_x,K33NLhat)
                calculateElement1!(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NLhat)
            end

        end

    end #END OF REDUCED INTEGRATION LOOP

    #unpack stored element stiffness data
    K11 = elStorage.K11
    K12 = elStorage.K12
    K13 = elStorage.K13
    K14 = elStorage.K14
    K15 = elStorage.K15
    K16 = elStorage.K16
    K22 = elStorage.K22
    K23 = elStorage.K23
    K24 = elStorage.K24
    K25 = elStorage.K25
    K26 = elStorage.K26
    K33 = elStorage.K33
    K34 = elStorage.K34
    K35 = elStorage.K35
    K36 = elStorage.K36
    K44 = elStorage.K44
    K45 = elStorage.K45
    K46 = elStorage.K46
    K55 = elStorage.K55
    K56 = elStorage.K56
    K66 = elStorage.K66

    if predef == "use"
        K21Pre = elStorage.K21NLpredef
        K12Pre = elStorage.K12NLpredef
        K31Pre = elStorage.K31NLpredef
        K13Pre = elStorage.K13NLpredef
        K22Pre = elStorage.K22NLpredef
        K23Pre = elStorage.K23NLpredef
        K33Pre = elStorage.K33NLpredef
    else
        K21Pre = zeros(size(K11))
        K12Pre = zeros(size(K11))
        K31Pre = zeros(size(K11))
        K13Pre = zeros(size(K11))
        K22Pre = zeros(size(K11))
        K23Pre = zeros(size(K11))
        K33Pre = zeros(size(K11))
    end

    if useDisp #modify stiffness matrices to account for nonlinear effects
        K21 = K12' + 2*K12NL' + K21Pre
        K12 = K12 + K12NL + K12Pre
        K31 = K13' + 2*K13NL' + K31Pre
        K13 = K13 + K13NL + K13Pre
        K22 = K22 + K22NL + K22Pre
        K23 = K23 + K23NL + K23Pre
        K33 = K33 + K33NL + K33Pre
    elseif useDisp == false && predef != "use"
        K21 = K12'
        K31 = K13'
    end

    if predef == "update"
        elStorage.K21NLpredef = collect(K21)
        elStorage.K12NLpredef = collect(K12)
        elStorage.K31NLpredef = collect(K31)
        elStorage.K13NLpredef = collect(K13)
        elStorage.K22NLpredef = collect(K22)
        elStorage.K23NLpredef = collect(K23)
        elStorage.K33NLpredef = collect(K33)
    end

    # Only used if (useDisp)
    K12hat  = K12
    K13hat  = K13

    #unpack stored element mass data
    M11 = elStorage.M11
    M15 = elStorage.M15
    M16 = elStorage.M16
    M22 = elStorage.M22
    M24 = elStorage.M24
    M33 = elStorage.M33
    M34 = elStorage.M34
    M44 = elStorage.M44
    M55 = elStorage.M55
    M56 = elStorage.M56
    M66 = elStorage.M66

    #unpack and scale stored element spin softening data
    S11 = elStorage.S11.*(O2^2 + O3^2)
    S12 = elStorage.S12.*O1*O2
    S13 = elStorage.S13.*O1*O3
    S15 = elStorage.S15.*(O2^2+O3^2)
    S16 = elStorage.S16.*(O2^2+O3^2)
    S22 = elStorage.S22.*(O1^2*0 + O3^2*0) #TODO: Why is this specific portion is causing increased nonphysical torque?
    S23 = elStorage.S23.*O2*O3
    S25 = elStorage.S25.*(O1*O2)
    S26 = elStorage.S26.*(O1*O2)
    S33 = elStorage.S33.*(O1^2+O2^2)
    S35 = elStorage.S35.*O1*O3
    S36 = elStorage.S36.*O1*O3
    S55 = elStorage.S55.*(O2^2+O3^2)
    S56 = elStorage.S56.*(O2^2+O3^2)
    S66 = elStorage.S66.*(O2^2+O3^2)
    S14 = elStorage.S14_1.*O1*O3 + elStorage.S14_2.*O1*O2
    S24 = elStorage.S24_1.*(O1^2+O3^2) + elStorage.S24_2.*O2*O3
    S34 = elStorage.S34_1.*(O1^2+O2^2) + elStorage.S34_2.*O2*O3
    S45 = elStorage.S45_1.*O1*O3 + elStorage.S45_2.*O1*O2
    S46 = elStorage.S46_1.*O1*O2 + elStorage.S46_2.*O1*O3
    S44 = elStorage.S44_1.*(O1^2+O3^2) + elStorage.S44_2.*(O1^2+O2^2) + elStorage.S44_3.*O2*O3

    #unpack and scale stored element Corilois data
    C12 = elStorage.C12.*O3
    C13 = elStorage.C13.*O2
    C23 = elStorage.C23.*O1
    C24 = elStorage.C24.*O1
    C25 = elStorage.C25.*O3
    C26 = elStorage.C26.*O3
    C34 = elStorage.C34.*O1
    C35 = elStorage.C35.*O2
    C36 = elStorage.C36.*O2
    C14 = elStorage.C14_1.*O2 + elStorage.C14_2.*O3
    C45 = elStorage.C45_1.*O3 + elStorage.C45_2.*O2
    C46 = elStorage.C46_1.*O2 + elStorage.C46_2.*O3

    #unpack and scale stored element Circulatory data
    H12 = 0.5*elStorage.C12.*O3dot
    H13 = 0.5*elStorage.C13.*O2dot
    H23 = 0.5*elStorage.C23.*O1dot
    H24 = 0.5*elStorage.C24.*O1dot
    H25 = 0.5*elStorage.C25.*O3dot
    H26 = 0.5*elStorage.C26.*O3dot
    H34 = 0.5*elStorage.C34.*O1dot
    H35 = 0.5*elStorage.C35.*O2dot
    H36 = 0.5*elStorage.C36.*O2dot
    H14 = 0.5*(elStorage.C14_1.*O2dot + elStorage.C14_2.*O3dot)
    H45 = 0.5*(elStorage.C45_1.*O3dot + elStorage.C45_2.*O2dot)
    H46 = 0.5*(elStorage.C46_1.*O2dot + elStorage.C46_2.*O3dot)


    #compile stiffness matrix without rotational effects
    Kenr = mapMatrixNonSym2(K11,K12,K13,K14,K15,K16,
    K21,K22,K23,K24,K25,K26,
    K31,K23',K33,K34,K35,K36,
    K13',K24',K34',K44,K45,K46,
    K15',K25',K35',K45',K55,K56,
    K16',K26',K36',K46',K56',K66)


    #add spin softening and circulatory effects to stiffness marix
    K11 = K11 .+ S11
    K21 = K21 .+ S12' .- H12'
    K12 = K12 .+ S12 .+ H12
    K31 = K31 .+ S13' .- H13'
    K13 = K13 .+ S13 .+ H13
    K41 = K14' .+ S14' .- H14'
    K14 = K14 .+ S14 .+ H14
    K15 = K15 .+ S15
    K16 = K16 .+ S16
    K22 = K22 .+ S22 .+ SS22
    K32 = K23' .+ S23' .- H23'
    K23 = K23 .+ S23 .+ H23
    K42 = K24'.+ S24' .- H24'
    K24 = K24 .+ S24 .+ H24
    K52 = K25'.+ S25' .- H25'
    K25 = K25 .+ S25 .+ H25
    K62 = K26' .+ S26' .- H26'
    K26 = K26 .+ S26 .+ H26
    K33 = K33 .+ S33 .+ SS33
    K43 = K34' .+ S34' .- H34'
    K34 = K34 .+ S34 .+ H34
    K53 = K35' .+ S35' .- H35'
    K35 = K35 .+ S35 .+ H35
    K63 = K36' .+ S36' .- H36'
    K36 = K36 .+ S36 .+ H36
    K44 = K44 .+ S44
    K54 = K45' .+ S45' .- H45'
    K45 = K45 .+ S45 .+ H45
    K64 = K46' .+ S46' .- H46'
    K46 = K46 .+ S46 .+ H46
    K55 = K55 .+ S55
    K56 = K56 .+ S56
    K66 = K66 .+ S66


    C43 = -C34'
    if (aeroElasticOn) #modify element matrices for aeroelastic effects
        K33 = K33 + K33Aero
        K34 = K34 + K34Aero
        C33 = C33 + C33Aero
        C34 = C34 + C34Aero
        K43 = K43 + K43Aero
        K44 = K44 + K44Aero
        C43 = C43 + C43Aero
        C44 = C44 + C44Aero
    end

    #---------------------------------------------
    zm=zeros(2,2)

    #compile stiffness matrix with rotational effects
    Ke = mapMatrixNonSym2(K11,K12,K13,K14,K15,K16,
    K21,K22,K23,K24,K25,K26,
    K31,K32,K33,K34,K35,K36,
    K41,K42,K43,K44,K45,K46,
    K15',K52,K53,K54,K55,K56,
    K16',K62,K63,K64,K56',K66)

    Kehat = zeros(12,12) # Declare type
    if useDisp && iterationType=="NR"
        #compile component of tangent matrix
        Kehat = mapMatrixNonSym2(zm,K12hat,K13hat,zm,zm,zm,
        zm,K22NLhat,K23NLhat,zm,zm,zm,
        zm,K23NLhat',K33NLhat,zm,zm,zm,
        zm,zm,zm,zm,zm,zm,
        zm,zm,zm,zm,zm,zm,
        zm,zm,zm,zm,zm,zm)
    end

    #compile Coriolis/damping matrix

    Ce = mapMatrixNonSym2(zm,C12,C13,C14,zm,zm,
    -C12',zm,C23,C24,C25,C26,
    -C13',-C23',C33,C34,C35,C36,
    -C14',-C24',C43,C44,C45,C46,
    zm,-C25',-C35',-C45',zm,zm,
    zm,-C26',-C36',-C46',zm,zm)

    #compile mass matrix
    Me = mapMatrixNonSym2(M11,zm,zm,zm,M15,M16,
    zm,M22,zm,M24,zm,zm,
    zm,zm,M33,M34,zm,zm,
    zm,M24',M34',M44,zm,zm,
    M15',zm,zm,zm,M55,M56,
    M16',zm,zm,zm,M56',M66)

    #account for rayleigh damping
    alpha = input.RayleighAlpha
    beta = input.RayleighBeta

    CeRayleigh = alpha.*Kenr + beta.*Me
    Ce = Ce + CeRayleigh

    #compile element force vector
    Fe = mapVector([F1;F2;F3;F4;F5;F6])

    # transform matrices for sweep
    # Note,a negative sweep angle, will sweep away from the direction of
    # positive rotation
    lambdaTran = lambda'

    Me = lambdaTran*Me*lambda
    Ce = lambdaTran*Ce*lambda

    Ke = lambdaTran*Ke*lambda
    if useDisp && iterationType=="NR"
        Kehat = lambdaTran*Kehat*lambda
    end

    Fe = lambdaTran*Fe

    ##concentrated mass
    #NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
    # if some ycm or zcm offset from the node was accounted for in concentrated mass terms
    ## Apply concentrated terms, including cross-coupling between concentrated mass and the other terms
    # TODO: should other cross-couplings be included here now since the off-diagonals can be included?

    concMassFlag = !isempty(findall(x->x!=0,concMass))
    concStiffFlag = !isempty(findall(x->x!=0,concStiff))
    concDampFlag = !isempty(findall(x->x!=0,concDamp))
    concLoadFlag = !isempty(findall(x->x!=0,concLoad))
    if concMassFlag
        #modify Me for concentrated mass
        Me[1:6,1:6] += concMass[:,1:6]
        Me[7:12,7:12] += concMass[:,7:12]
    end

    if concMassFlag || concDampFlag
        Ce[1:6,1:6] += concDamp[:,1:6]
        Ce[7:12,7:12] += concDamp[:,7:12]

        #modify Ce for concentrated mass
        Ce[1,2] = Ce[1,2] - 2*concMass[1,1]*omega_z
        Ce[2,1] = Ce[2,1] + 2*concMass[1,1]*omega_z
        Ce[1,3] = Ce[1,3] + 2*concMass[1,1]*omega_y
        Ce[3,1] = Ce[3,1] - 2*concMass[1,1]*omega_y
        Ce[2,3] = Ce[2,3] - 2*concMass[1,1]*omega_x
        Ce[3,2] = Ce[3,2] + 2*concMass[1,1]*omega_x
        Ce[7,8] = Ce[7,8] - 2*concMass[1,7]*omega_z
        Ce[8,7] = Ce[8,7] + 2*concMass[1,7]*omega_z
        Ce[7,9] = Ce[7,9] + 2*concMass[1,7]*omega_y
        Ce[9,7] = Ce[9,7] - 2*concMass[1,7]*omega_y
        Ce[8,9] = Ce[8,9] - 2*concMass[1,7]*omega_x
        Ce[9,8] = Ce[9,8] + 2*concMass[1,7]*omega_x
    end

    if concMassFlag || concStiffFlag
        #modify Ke for concentrated mass
        Ke[1:6,1:6] += concStiff[:,1:6]
        Ke[7:12,7:12] += concStiff[:,7:12]
        Ke[1,1] -= concMass[1,1]*(omega_y^2 + omega_z^2)
        Ke[1,2] += concMass[1,1]*omega_x*omega_y - concMass[1,1]*omegaDot_z
        Ke[2,1] += concMass[1,1]*omega_x*omega_y + concMass[1,1]*omegaDot_z
        Ke[1,3] += concMass[1,1]*omega_x*omega_z + concMass[1,1]*omegaDot_y
        Ke[3,1] += concMass[1,1]*omega_x*omega_z - concMass[1,1]*omegaDot_y
        Ke[2,3] += concMass[1,1]*omega_y*omega_z - concMass[1,1]*omegaDot_x
        Ke[3,2] += concMass[1,1]*omega_y*omega_z + concMass[1,1]*omegaDot_x
        Ke[2,2] -= concMass[1,1]*(omega_x^2 + omega_z^2)
        Ke[3,3] -= concMass[1,1]*(omega_x^2 + omega_y^2)
        Ke[7,7] -= concMass[1,7]*(omega_y^2 + omega_z^2)
        Ke[7,8] += concMass[1,7]*omega_x*omega_y - concMass[1,7]*omegaDot_z
        Ke[8,7] += concMass[1,7]*omega_x*omega_y + concMass[1,7]*omegaDot_z
        Ke[7,9] += concMass[1,7]*omega_x*omega_z + concMass[1,7]*omegaDot_y
        Ke[9,7] += concMass[1,7]*omega_x*omega_z - concMass[1,7]*omegaDot_y
        Ke[8,9] += concMass[1,7]*omega_y*omega_z - concMass[1,7]*omegaDot_x
        Ke[9,8] += concMass[1,7]*omega_y*omega_z + concMass[1,7]*omegaDot_x
        Ke[8,8] -= concMass[1,7]*(omega_x^2 + omega_z^2)
        Ke[9,9] -= concMass[1,7]*(omega_x^2 + omega_y^2)
    end

    if concMassFlag || concLoadFlag
        Fe[1:6] += concLoad[:,1]
        Fe[7:12] += concLoad[:,2]

        #modify Fe for  concentrated load
        Fe[1] += concMass[1,1]*(x[1]*(omega_y^2 + omega_z^2)-omega_x*omega_y*y[1] - omega_x*omega_z*z[1]) + concMass[1,1]*(y[1]*omegaDot_z-z[1]*omegaDot_y)  -  concMass[1,1]*a_x
        Fe[2] += concMass[1,1]*(y[1]*(omega_x^2 + omega_z^2)-omega_y*omega_z*z[1] - omega_y*omega_x*x[1]) + concMass[1,1]*(z[1]*omegaDot_x-x[1]*omegaDot_z)  -  concMass[1,1]*a_y
        Fe[3] += concMass[1,1]*(z[1]*(omega_x^2 + omega_y^2)-omega_z*omega_x*x[1] - omega_z*omega_y*y[1]) + concMass[1,1]*(x[1]*omegaDot_y-y[1]*omegaDot_x)  -  concMass[1,1]*a_z
        Fe[7] += concMass[1,7]*(x[2]*(omega_y^2 + omega_z^2)-omega_x*omega_y*y[2] - omega_x*omega_z*z[2]) + concMass[1,7]*(y[2]*omegaDot_z-z[2]*omegaDot_y)  -  concMass[1,7]*a_x
        Fe[8] += concMass[1,7]*(y[2]*(omega_x^2 + omega_z^2)-omega_y*omega_z*z[2] - omega_y*omega_x*x[2]) + concMass[1,7]*(z[2]*omegaDot_x-x[2]*omegaDot_z)  -  concMass[1,7]*a_y
        Fe[9] += concMass[1,7]*(z[2]*(omega_x^2 + omega_y^2)-omega_z*omega_x*x[2] - omega_z*omega_y*y[2]) + concMass[1,7]*(x[2]*omegaDot_y-y[2]*omegaDot_x)  -  concMass[1,7]*a_z
    end


    ##
    # Declare Types
    Fhate = zeros(12)
    FhatLessConc = zeros(12)
    if analysisType=="TD" #calculate effective stiffness matrix and force vector for Dean integrator

        a1 = timeInt.a1
        a2 = timeInt.a2
        a3 = timeInt.a3
        a4 = timeInt.a4

        xn=disp[1:12]
        xnm1=dispm1[1:12]
        A = 2.0.*xn - xnm1
        B = -a1.*xnm1 - a2.*xn
        D = a3.*xnm1

        Khate = Ke*a1 + a3.*Ce + Me
        Fhate = Fe*a4 + Me*(A) + Ke*(B) + Ce*(D)

        FhatLessConc =   Fhate - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]].*a4

        #........................................................

        #..........................................................

        Ke = copy(Khate)
        Fe = copy(Fhate)
    end

    if analysisType=="TNB" #calculate effective stiffness matrix and load vector for Newmark-Beta integrator
        #     a1 = timeInt.a1
        #     a2 = timeInt.a2
        a3 = timeInt.a3
        a4 = timeInt.a4
        a5 = timeInt.a5
        a6 = timeInt.a6
        a7 = timeInt.a7
        a8 = timeInt.a8

        u=copy(disp)
        udot=copy(dispdot)
        uddot=copy(dispddot)
        if (iterationType=="NR")    #considerations if newton raphson iteration is used
            if (input.firstIteration)
                A = a3*u + a4*udot + a5*uddot
                B = a6*u + a7*udot + a8*uddot
                Fhate = Fe + Me*(A) + Ce*(B) - Ke*u
            else
                Fhate = Fe  - Me*uddot - Ce*udot - (Ke)*u
            end
        elseif iterationType=="DI"||iterationType=="LINEAR"   #considerations if direct iteration is used or linear analysis
            A = a3*u + a4*udot + a5*uddot
            B = a6*u + a7*udot + a8*uddot
            Fhate = Fe + Me*(A) + Ce*(B)
        end

        Khate = Ke + a3.*Me + a6.*Ce
        if iterationType=="NR" #considerations if newton raphson iteration is used
            Khate = Kehat + Khate
        end

        FhatLessConc =   Fhate - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]]

        #........................................................

        Ke = copy(Khate)
        Fe = copy(Fhate)

    end

    if (analysisType=="M")
        FhatLessConc =   Fe - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]]

        if (iterationType=="DI")
            Fe = Fe*input.loadStep
        end
    end

    if ((analysisType=="M" || analysisType=="S") && iterationType=="NR") #considerations for newton-raphson iteration
        Fe = Fe*input.loadStep - Ke*disp_iter
        Ke = Ke .+ Kehat
    end

    ###----- assign output block ----------------

    if !(analysisType=="M"||analysisType=="ROM")
        Me = zeros(size(Me))
        Ce = zeros(size(Ce))
    end

    return ElOutput(FhatLessConc,Ke,Fe,Me,Ce)
end

"""
    calculateTimoshenkoElementStrain(elementOrder,nlOn,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,disp)

Internal, calculates element strain for a Timoshenko element

#Outputs
* `ElStrain`: See ?ElStrain

"""
function calculateTimoshenkoElementStrain(elementOrder,nlOn,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,disp)

    numGP = 4   #number of gauss points for full integration
    #calculate quad points
    xi,_ = getGP(numGP)

    # p_disp_x = zeros(numGP,6)

    #Initialize element sub matrices and sub vectors

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg_d = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg_d.*pi/180.0)

    dispLocal = lambda*disp'

    uNode = [dispLocal[1] dispLocal[7]]
    vNode = [dispLocal[2] dispLocal[8]]
    wNode = [dispLocal[3] dispLocal[9]]
    theta_xNode = [dispLocal[4]  dispLocal[10]]
    theta_yNode = [dispLocal[5]  dispLocal[11]]
    theta_zNode = [dispLocal[6]  dispLocal[12]]

    #Integration loop
    epsilon_x = zeros(numGP)
    epsilon_y = zeros(numGP)
    epsilon_z = zeros(numGP)
    kappa_x = zeros(numGP)
    kappa_y = zeros(numGP)
    kappa_z = zeros(numGP)
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,_ = calculateShapeFunctions(elementOrder,xi[i],xloc)
        #N1 = N
        #N2 = N
        #N3 = N
        #N4 = N
        N5 = copy(N)
        N6 = copy(N)
        p_N1_x = copy(p_N_x)
        p_N2_x = copy(p_N_x)
        p_N3_x = copy(p_N_x)
        p_N4_x = copy(p_N_x)
        p_N5_x = copy(p_N_x)
        p_N6_x = copy(p_N_x)

        #calculate displacement derivatives at quad point i
        uprime = interpolateVal(uNode,p_N1_x)
        vprime = interpolateVal(vNode,p_N2_x)
        wprime = interpolateVal(wNode,p_N3_x)
        theta_x_prime = interpolateVal(theta_xNode,p_N4_x)
        theta_y_prime = interpolateVal(theta_yNode,p_N5_x)
        theta_y_gp = interpolateVal(theta_yNode,N5)
        theta_z_prime = interpolateVal(theta_zNode,p_N6_x)
        theta_z_gp = interpolateVal(theta_zNode,N6)
        # p_disp_x[i,:] = [uprime, vprime, wprime, theta_x_prime, theta_y_prime, theta_z_prime]

        if nlOn
            epsilon_x[i] = uprime + 0.5*(wprime^2 + vprime^2) #e_xx_0
        else
            epsilon_x[i] = uprime
        end
        epsilon_y[i] = -theta_z_gp + vprime #g_xy_0
        epsilon_z[i] =  theta_y_gp + wprime #g_xz_0
        kappa_x[i] =  theta_x_prime #g_xz_y
        kappa_y[i] = theta_y_prime #e_xx_z
        kappa_z[i] = -theta_z_prime #e_xx_y
    end #END OF INTEGRATION LOOP

    return ElStrain(epsilon_x,epsilon_y,epsilon_z,kappa_x,kappa_y,kappa_z)
end

"""

    calculateTimoshenkoElementNLSS(elinput)

Performs selective nonlinear element calculations. Only stiffness matrix contributions are evaluate. No other calculations are performed to facilitate efficiency gains.

#Input
*  `elinput`:     object containing element input

#Output
*  `eloutput`:    object containing element data
"""
function calculateTimoshenkoElementNLSS(input)

    ###-------- assign input block ----------------
    elementOrder   = input.elementOrder
    x              = input.x
    xloc           = input.xloc
    disp           = input.disp
    sectionProps   = input.sectionProps
    sweepAngle     = input.sweepAngle
    coneAngle      = input.coneAngle
    rollAngle      = input.rollAngle

    useDisp        = input.useDisp
    preStress      = input.preStress
    iterationType  = input.iterationType

    ###--------------------------------------------
    if input.analysisType == "M" #TODO: why are we doing this if the analysis type is hard coded above to be M and required below?
        disp_iter=disp
    end

    numGP = 1 #used reduced integration for nonlinear terms

    #calculate quad points
    xi,weight = getGP(numGP)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    K22NL = zeros(numNodesPerEl,numNodesPerEl)
    K33NL = zero(K22NL)
    K12NL = zero(K22NL)
    K13NL = zero(K22NL)
    K23NL = zero(K22NL)
    #     SS33 = SS22

    #Convert frequencies from Hz to radians
    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    dispLocal = lambda*disp_iter[1:12]

    uNode = [dispLocal[1] dispLocal[7]]
    vNode = [dispLocal[2] dispLocal[8]]
    wNode = [dispLocal[3] dispLocal[9]]

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        p_N1_x = p_N_x
        p_N2_x = p_N_x
        p_N3_x = p_N_x

        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        if useDisp || preStress
            EA   = interpolateVal(sectionProps.EA,N)
            uprime = interpolateVal(uNode,p_N1_x)
            vprime = interpolateVal(vNode,p_N2_x)
            wprime = interpolateVal(wNode,p_N3_x)
        end

        #nonlinear element calculations
        if useDisp
            calculateElement1!(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12NL)
            calculateElement1!(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13NL)
            calculateElement1!(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22NL)
            calculateElement1!(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NL)
            calculateElement1!(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33NL)

            #K12NLhat = K12
            #K13NLhat = K13
            #nonlinear element tangent matrix component calculations
            # T_ij = K_ij + Khat_ij
            if iterationType == "NR"
                calculateElement1!(EA*(uprime + vprime^2 + 0.5*wprime^2),integrationFactor,p_N2_x,p_N2_x,K22NLhat)
                calculateElement1!(EA*(uprime + wprime^2 + 0.5*vprime^2),integrationFactor,p_N3_x,p_N3_x,K33NLhat)
                calculateElement1!(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NLhat)
            end

        end


    end #END OF INTEGRATION LOOP

    ###---------------------------------------------

    Ktemp = zeros(numNodesPerEl*6,numNodesPerEl*6)
    Ktemp[1:2,3:4] = K12NL
    Ktemp[1:2,5:6] = K13NL
    Ktemp[3:4,1:2] = 2*K12NL
    Ktemp[5:6,1:2] = 2*K13NL
    Ktemp[3:4,3:4] = K22NL
    Ktemp[5:6,5:6] = K33NL
    Ktemp[3:4,5:6] = K23NL
    Ktemp[5:6,3:4] = K23NL

    Ke = mapMatrixNonSym(Ktemp)

    # transform matrices for sweep (currently hardcoded to 0 sweep)
    # Note,a negative theta3, will sweep away from the direction of
    # rotation
    lambdaTran = lambda'
    Ke = SparseArrays.sparse(Ke)
    lambda = SparseArrays.sparse(lambda)
    lambdaTran = SparseArrays.sparse(lambdaTran)
    Ke = lambdaTran*Ke*lambda

    if iterationType == "NR"
        error("calcTimoElNLSS needs some mods to be used with newton raphson")
    end
    #----- assign output block ----------------
    Ke = collect(Ke)

    ###------------------------------------------
    return Ke
end
