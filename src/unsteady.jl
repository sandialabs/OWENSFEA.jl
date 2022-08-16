
"""
    initialElementCalculations(feamodel,el,mesh)

performs initial element calculation for use later in analysis for efficiency gains.

# Inputs
* `feamodel::FEAmodel`: see ?Feamodel
* `el::El`: see ?El
* `mesh::Mesh`: see ?Mesh

# Outputs:
* `elStorage:ElStorage`: see ?ElStorage

"""
function initialElementCalculations(feamodel,el,mesh)

    #initial element calculation
    numNodesPerEl = 2

    elStorage = Array{ElStorage, 1}(undef, mesh.numEl)
    elx = zeros(numNodesPerEl)
    ely = zeros(numNodesPerEl)
    elz = zeros(numNodesPerEl)

    elementOrder = feamodel.elementOrder #assign for element i
    modalFlag = true
    aeroSweepAngle = 0.0
    for i=1:mesh.numEl
        #Calculate Ke and Fe for element i
        xloc = [0.0 el.elLen[i]]
        sectionProps = el.props[i]
        sweepAngle = el.psi[i]
        coneAngle = el.theta[i]
        rollAngle = el.roll[i]

        for j=1:numNodesPerEl
            #get element cooridnates
            elx[j] = mesh.x[Int(mesh.conn[i,j])]
            ely[j] = mesh.y[Int(mesh.conn[i,j])]
            elz[j] = mesh.z[Int(mesh.conn[i,j])]
        end

        #get concentrated terms associated with element
        massConc,_,_,_,_,_ = ConcMassAssociatedWithElement(mesh.conn[i,:],feamodel.joint,feamodel.nodalTerms.concMass,feamodel.nodalTerms.concStiff,feamodel.nodalTerms.concLoad)

        concMassFlag = !isempty(findall(x->x!=0,massConc))

        Omega = 0.0

        elStorage[i] = calculateTimoshenkoElementInitialRun(elementOrder,modalFlag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,elx,ely,elz,concMassFlag,massConc,Omega) #initial element calculations for storage

    end
    return elStorage
end

"""
    structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega,OmegaDot,time,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)

performs unsteady structural dynamics analysis

# Inputs
* `feamodel::`: object containing feamodel data
* `mesh::`: object containing mesh data
* `el::`: object containing element data
* `dispData::`: object containing displacement data
* `Omega::`: rotor speed (Hz)
* `OmegaDot::`: rotor acceleratin (Hz)
* `time::`: current simulation time
* `delta_t::`: time step size
* `elStorage::`: object containing stored element data
* `Fexternal::`: vector containing external force values
* `Fdof::`: vector containing global DOF numbering associated with external force values
* `CN2H::`: transformation matrix from inertial frame to hub frame
* `rbData::`: vector containing rigid body displacement, velocity, and acceleration

# Outputs:
* `elStrain::ElStrain`: see ?ElStrain strain for element at end of time step
* `dispOut::DispOut`: see ?DispOut displacement data at end of time step
* `FReaction_sp1::`: vector containing reaction force at turbine base at end of time step
"""
function  structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega,OmegaDot,time,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData;predef=nothing)

    #-------- get feamodel information -----------
    conn = Int.(mesh.conn) #TODO: make the generator output ints
    iterationType = feamodel.nlParams.iterationType #"NR" "DI" "LINEAR"
    analysisType = feamodel.analysisType

    #initialize displacements, tolerance, uNorm, iteration count for nonlinear iteration
    unorm = 1e6
    tol = feamodel.nlParams.tolerance
    maxIterations = feamodel.nlParams.maxIterations
    iterationCount = 0

    # Initialize data which will be iterated on
    elementOrder = feamodel.elementOrder
    numNodesPerEl = elementOrder + 1
    numDOFPerNode = 6
    totalNumDOF = mesh.numNodes * numDOFPerNode
    eldisp = zeros(numNodesPerEl*numDOFPerNode)
    eldisp_sm1 = zeros(numNodesPerEl*numDOFPerNode)
    Kg1 = zeros(totalNumDOF,totalNumDOF) #initialize global stiffness and force vector
    Fg1 = zeros(totalNumDOF)
    eldispdot = zero(eldisp)
    eldispddot = zero(eldisp)
    eldispiter = zero(eldisp)
    timeInt = TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) # initialize, is filled in later

    disp_s = copy(dispData.displ_s)
    dispdot_s = copy(dispData.displdot_s)
    dispddot_s = copy(dispData.displddot_s)
    disp_sm1 = copy(dispData.displ_sm1)

    displddot_im1 = copy(dispddot_s)
    displdot_im1 = copy(dispdot_s)
    displ_im1 = copy(disp_s)

    while (unorm>tol && iterationCount < maxIterations) #iteration loop
        if iterationCount>maxIterations-5
            println("$unorm $iterationCount")
        end

        Kg = zero(Kg1)
        Fg = zero(Fg1)
        nodalTerms,timeInt = TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,
        dispData,Omega,elStorage;Kg,Fg,eldisp_sm1,eldispdot,eldispddot,eldispiter,rbData,CN2H,delta_t,
        OmegaDot,displ_im1,displdot_im1,displddot_im1,iterationCount,predef)

        #apply general 6x6  mass, damping, and stiffness matrices to nodes
        # Except only stiffness is used here... #TODO: is this correct?  Shouldn't the cross coupling between say mass and force (from acceleration) be included with the cross terms?
        # The way it's coded, if a 6x6 matrix is used, the diagonals are used in the timoshenko dynamic elements, and the cross terms (without the diagonal) are applied here.
        Kg,_,_ = applyGeneralConcentratedTerms(Kg,Kg,Kg,nodalTerms.concStiffGen,nodalTerms.concMassGen,nodalTerms.concDampGen)

        #Apply external loads to structure
        for i=1:length(Fexternal)
            if analysisType=="TD"
                Fg[Fdof[i]] = Fg[Fdof[i]] + Fexternal[i]*delta_t^2
            end
            if analysisType=="TNB"
                Fg[Fdof[i]] = Fg[Fdof[i]] + Fexternal[i]
            end
        end

        #------ apply constraints on system -----------------------------------
        Kg = applyConstraints(Kg,feamodel.jointTransform)
        Fg = applyConstraintsVec(Fg,feamodel.jointTransform)

        #Apply BCs to global system
        KgTotal,FgTotal = applyBC(Kg,Fg,feamodel.BC,numDOFPerNode)
        solution = KgTotal\FgTotal  #solve for displacements

        solution = feamodel.jointTransform*solution #transform to full dof listing

        if feamodel.nlOn  #calculate norm between current iteration and last iteration
            if iterationType=="NR"
                unorm = calcUnorm(displ_im1+solution,displ_im1)
            else
                unorm = calcUnorm(solution,displ_im1)
            end
        else
            unorm = 0.0
        end

        if iterationType=="NR"
            #if newton raphson update u, udot, uddot at each iteration
            displ_im1 = displ_im1 + solution
            cap_delta_displ = displ_im1 - dispData.displ_s
            displddot_im1 = timeInt.a3*(cap_delta_displ) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s
            displdot_im1  = -timeInt.a7*dispData.displdot_s -timeInt.a8*dispData.displddot_s + timeInt.a6*(cap_delta_displ)
        elseif (iterationType=="DI"||iterationType=="LINEAR")
            displ_im1 = solution
        else
            error("iteration type not supported, choose another")
        end

        iterationCount = iterationCount + 1
    end #While

    #Calculate reaction at turbine base (hardwired to node number 1)
    reactionNodeNumber = feamodel.platformTurbineConnectionNodeNumber

    FReaction = calculateReactionForceAtNode(reactionNodeNumber,feamodel,mesh,el,elStorage,timeInt,dispData,displ_im1,rbData,Omega,OmegaDot,CN2H)
    #Calculate strain
    elStrain = calculateStrainForElements(mesh.numEl,numNodesPerEl,numDOFPerNode,conn,feamodel.elementOrder,el,displ_im1,feamodel.nlOn)
    if (iterationCount>=maxIterations)
        @warn "Maximum iterations exceeded."
    end

    FReaction_sp1 = FReaction
    displ_sp1 = displ_im1

    # Specific to TNB, but must be declared
    displddot_sp1 = timeInt.a3*(displ_sp1-dispData.displ_s) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s #store acceleration vector in dispOut
    displdot_sp1 = dispData.displdot_s + timeInt.a2*dispData.displddot_s + timeInt.a1*displddot_sp1                    #store velocity vector in dispOut

    dispOut = DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1)

    return elStrain,dispOut,FReaction_sp1
end

"""
Internal, this function transforms a matrix by the transformation matrix to enforce joint constraints
"""
function applyConstraints(Kg,transMatrix)
    return transMatrix'*(Kg*transMatrix)
end

"""
Internal, this function transforms a vector by the transformation matrix to enforce joint constraints
"""
function applyConstraintsVec(Fg,transMatrix)
    return transMatrix'*Fg
end

"""
This function calculates a relative norm between two vectors: unew and uold
"""
function calcUnorm(unew,uold)
    return LinearAlgebra.norm(unew-uold)/LinearAlgebra.norm(unew)
end

"""
Internal, function to form total stifness matrix and transform to desired DOF mapping
"""
function mapMatrixNonSym2(K11,K12,K13,K14,K15,K16,K21,K22,K23,K24,K25,K26,K31,K32,K33,K34,K35,K36,K41,K42,K43,K44,K45,K46,K51,K52,K53,K54,K55,K56,K61,K62,K63,K64,K65,K66)

    T = [1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

    Ktemp = zeros(12,12)

    Ktemp[1:2,1:2] = K11
    Ktemp[1:2,3:4] = K12
    Ktemp[1:2,5:6] = K13
    Ktemp[1:2,7:8] = K14
    Ktemp[1:2,9:10] = K15
    Ktemp[1:2,11:12] = K16

    Ktemp[3:4,1:2] = K21
    Ktemp[3:4,3:4] = K22
    Ktemp[3:4,5:6] = K23
    Ktemp[3:4,7:8] = K24
    Ktemp[3:4,9:10] = K25
    Ktemp[3:4,11:12] = K26

    Ktemp[5:6,1:2] = K31
    Ktemp[5:6,3:4] = K32
    Ktemp[5:6,5:6] = K33
    Ktemp[5:6,7:8] = K34
    Ktemp[5:6,9:10] = K35
    Ktemp[5:6,11:12] = K36

    Ktemp[7:8,1:2] = K41
    Ktemp[7:8,3:4] = K42
    Ktemp[7:8,5:6] = K43
    Ktemp[7:8,7:8] = K44
    Ktemp[7:8,9:10] = K45
    Ktemp[7:8,11:12] = K46

    Ktemp[9:10,1:2] = K51
    Ktemp[9:10,3:4] = K52
    Ktemp[9:10,5:6] = K53
    Ktemp[9:10,7:8] = K54
    Ktemp[9:10,9:10] = K55
    Ktemp[9:10,11:12] = K56

    Ktemp[11:12,1:2] = K61
    Ktemp[11:12,3:4] = K62
    Ktemp[11:12,5:6] = K63
    Ktemp[11:12,7:8] = K64
    Ktemp[11:12,9:10] = K65
    Ktemp[11:12,11:12] = K66

    #map to FEA numbering
    # T2 = SparseArrays.sparse(T)
    Kel = T'*Ktemp*T

    #declare map
    # map = [1, 7, 2, 8, 3, 9,...
    #       4, 10, 5, 11, 6, 12];
    #
    # #map to FEA numbering
    # for i=1:a
    #     I=map[i];
    #     for j=1:a
    #         J=map(j);
    #         Kel(I,J) = Ktemp(i,j);
    #     end
    # end

    return Kel

end

"""
Internal, function to form total stifness matrix and transform to desired DOF mapping
"""
function mapMatrixNonSym(Ktemp)

    T = [1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

    #map to FEA numbering
    T2 = SparseArrays.sparse(T)
    Kel = T2'*Ktemp*T2

    #declare map
    # map = [1, 7, 2, 8, 3, 9,...
    #       4, 10, 5, 11, 6, 12];
    #
    # #map to FEA numbering
    # for i=1:a
    #     I=map[i];
    #     for j=1:a
    #         J=map(j);
    #         Kel(I,J) = Ktemp(i,j);
    #     end
    # end

    return Kel

end

"""
    elementPostProcess(elementNumber,feamodel,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)

Internal, calculates the reaction force associated with an element.

#Input
* `elementNumber`  node number joint constraints are desired at
* `feamodel`       object containing feamodel data
* `mesh`           object containing mesh data
* `elStorage`      object containing stored element data
* `el`             object containing element data
* `timeInt`        object containing time integration parameters
* `dispData`       object containing displacement data
* `displ_iter`     converged displacement solution
* `rbData`         vector containing rigid body displacement, velocity, and acceleration
* `Omega`          rotor speed (Hz)
* `OmegaDot`       rotor acceleratin (Hz)
* `CN2H`           transformation matrix from inertial frame to hub frame

#Output
* `Fpp`            vector containing reaction force vector associated with element
"""
function elementPostProcess(elementNumber,feamodel,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
    elementOrder = feamodel.elementOrder
    numNodesPerEl = elementOrder + 1
    numDOFPerNode = 6
    eldisp = zeros(numNodesPerEl*numDOFPerNode)

    Fpp = TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,dispData,Omega,elStorage;
        rbData,CN2H,OmegaDot,displ_im1=displ_iter,timeInt,postprocess=true,elementNumber)

    return Fpp
end
