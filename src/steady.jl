"""

staticAnalysis(feamodel,mesh,el,displ,Omega,OmegaStart,elStorage;
    reactionNodeNumber=1, OmegaDot=0.0, Fdof=[1], Fexternal=[0.0])

This function performs a static analysis and returns displacement
values and a flag denoting successful/unsuccessful analysis

#Inputs
* `feamodel`:          object containing feamodel information
* `mesh`:           object containing mesh information
* `el`:             object containing element information
* `displ`:          displacement vector for use in pre-stressed analysis
* `Omega`:          rotor speed (Hz)
* `OmegaStart`:     rotor speed (Hz) from previous analysis if stepping through various rotor speeds, may be useful in load stepping
* `elStorage`:      previously calculated element system matrices
* `reactionNodeNumber::Int`: optional, node at which to calculate reaction force
* `OmegaDot::Float`: Steady State Rotational Acceleration
* `Fdof::Array{<:Int}`: Global Dofs where Fexternal is acting, where max dof = nelem*ndof
* `Fexternal{<:Float}`: Forces or moments associated with the Fdofs specified
#Outputs
* `displ`:                    vector of displacemetns
* `staticAnalysisSuccessful`: boolean flag denoting successful static analysis
"""
function staticAnalysis(feamodel,mesh,el,displ,Omega,OmegaStart,elStorage;
    reactionNodeNumber=1, OmegaDot=0.0, Fdof=[1], Fexternal=[0.0])

    feamodel.analysisType = "S" #Force type to align with the static/steady call

    # x = mesh.x
    # y = mesh.y
    # z = mesh.z
    conn = mesh.conn
    # omegaVec       = zeros(3)
    # omegaDotVec    = zeros(3)
    # accelVec       = zeros(3)
    # CN2H           = 1.0*LinearAlgebra.I(3)
    # dispm1 = zeros(12) #declare type
    dispOld = copy(displ) #initialize scope
    staticAnalysisSuccessfulForLoadStep = false #initialize scope
    countedNodes = [] #TODO:??

    elementOrder = feamodel.elementOrder #extract element order from feamodel
    numNodesPerEl = elementOrder + 1 #do initialization
    numDOFPerNode = 6
    totalNumDOF = mesh.numNodes * numDOFPerNode

    # nodalTerms = feamodel.nodalTerms #extract concentrated nodal terms from feamodel
    # nodalTermsCopy = deepcopy(nodalTerms)   #extract extra copy of concentrated nodal terms from feamodel

    #load stepping paramters
    loadStepCount = 1
    displPrev = copy(displ) #copy of initial displacement vector
    staticAnalysisComplete = false #initialize to false
    nlParams = feamodel.nlParams

    if nlParams.adaptiveLoadSteppingFlag
        loadStepPrev = 0.0
        loadStep = 1.0
    else
        loadStepPrev = 0.0
        loadStep = nlParams.prescribedLoadStep[1]
        println("Prescribed load step: $loadStep")
    end

    maxNumLoadSteps = nlParams.maxNumLoadSteps
    MAXIT = nlParams.maxIterations
    tolerance = nlParams.tolerance
    # timeInt = TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) #not used
    eldisp = zeros(numNodesPerEl*numDOFPerNode)
    Kg1 = zeros(totalNumDOF,totalNumDOF)   #initialize global stiffness matrix
    Fg1 = zeros(totalNumDOF)             #initialize global force vector
    #.........................................................................
    while !staticAnalysisComplete && loadStepCount<maxNumLoadSteps
        # staticAnalysisSuccessful = false #initialize staticAnalysisSuccessful flag
        uNorm = 1.0 #initialize norm for convergence check

        iterationCount = 0 #initialize iteration count

        while uNorm > tolerance && iterationCount < MAXIT #iteration loop (convergence tolerance of 1.0e-6)
            Kg = zero(Kg1)
            Fg = zero(Fg1)
            TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,displ,Omega,elStorage;
                Kg,Fg,iterationCount,dispOld,loadStepPrev,loadStep,OmegaDot,countedNodes)

            # Fexternal, Fdof = externalForcingStatic()  #TODO: get arbitrary external loads from externalForcingStatic() function
            for i=1:length(Fdof)
                Fg[Fdof[i]] =  Fg[Fdof[i]] + Fexternal[i]*loadStep #modify assembled global load vector for external loads
            end

            #----------------------------------------------------------------------

            # #apply general 6x6  mass, damping, and stiffness matrices to nodes
            # Kg,_,_ = applyGeneralConcentratedTerms(Kg,Kg,Kg,feamodel.nodalTerms.concStiffGen,feamodel.nodalTerms.concMassGen,feamodel.nodalTerms.concDampGen)

            #APPLY BOUNDARY CONDITIONS
            Kg = applyConstraints(Kg,feamodel.jointTransform) #modify global stiffness matrix for joint constraints using joint transform
            Fg = applyConstraintsVec(Fg,feamodel.jointTransform) #modify global force vector for joint constraints using joint transform

            if feamodel.BC.numpBC==0
                @warn "No boundary conditions detected. Fully fixing DOFs at Node 1 to faciliate static solve."
                feamodel.BC.numpBC = 6
                feamodel.BC.pBC = [1 1 0; 1 2 0; 1 3 0;1 4 0;1 5 0;1 6 0];
            end

            Kg,Fg = applyBC(Kg,Fg,feamodel.BC,numDOFPerNode)  #apply boundary conditions to global stiffness matrix and force vector
            dispOld = copy(displ)  #assign displacement vector from previous iteration

            if nlParams.iterationType == "NR"  #system solve, norm calculation for newton-raphson iteration
                delta_displ = Kg\Fg
                delta_displ = feamodel.jointTransform*delta_displ
                displ = displ + delta_displ
                uNorm = calcUnorm(displ-delta_displ,displ)
            elseif nlParams.iterationType == "DI" #system solve, norm calculation for direct iteration
                displ_last = copy(displ)
                displ = Kg\Fg
                displ = feamodel.jointTransform*displ
                uNorm = calcUnorm(displ_last,displ)
                gamm = 0.5
                displ = (1-gamm)*displ + gamm*displ_last
            else                                        #system solve for linear case
                displ = Kg\Fg
                displ = feamodel.jointTransform*displ
                uNorm = 0
            end
            iterationCount = iterationCount +1         #increment iteration count
        end

        loadStepCount = loadStepCount + 1 #increment load step count

        #update load step whether adaptive or prescribed
        loadStep,loadStepPrev,displ,displPrev,staticAnalysisSuccessfulForLoadStep,staticAnalysisComplete = updateLoadStep(iterationCount,nlParams,loadStep,loadStepPrev,loadStepCount,displPrev,displ)

    end
    staticAnalysisSuccessful = staticAnalysisSuccessfulForLoadStep
    # t_static = toc
    # println("Elapsed time for static analysis(s):")
    # println(t_static)

    elStrain = calculateStrainForElements(mesh.numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ,feamodel.nlOn)  
    #feamodel.platformTurbineConnectionNodeNumber #TODO: multiple points?  the whole mesh?
    timeInt = nothing
    dispData = copy(displ)
    rbData = zeros(9)
    CN2H = 1.0*LinearAlgebra.I(3)
    #Calculate reaction at turbine base (hardwired to node number 1)
    FReaction = zeros(mesh.numEl*6)
    for reactionNodeNumber = 1:mesh.numEl
        try
            countedNodes = [] #TODO:??
            FReaction[(reactionNodeNumber-1)*6+1:reactionNodeNumber*6] = calculateReactionForceAtNode(reactionNodeNumber,feamodel,mesh,el,elStorage,timeInt,dispData,displ,rbData,Omega,OmegaDot,CN2H,countedNodes;single_element_reaction=true)
        catch
            # This is where a joint is println(reactionNodeNumber)
        end
    end
    
    return displ,elStrain,staticAnalysisSuccessful,FReaction

end

"""

    updateLoadStep(iterationCount,loadStepParams,loadStep,loadStepPrev,loadStepCount,displCopy,displ)

Updates the load stepping parameter whether through means
of adaptive loadstepping or a specified load step profile.

#Input
* `iterationCount`      number of iterations for current load step
* `loadStepParams`      struct containing load step parameters
* `loadStep`            load step value for current load step
* `loadStepPrev`        load step value for previous load st ep
* `loadStepCount`       number of load steps performed up to this point
* `displPrev`           converged displacement vector form previous load step
* `displ`               displacement vector at current load step

#Output
* `loadStep`            new load step value
* `loadStepPrev`        load step value for previous load step
* `displ`               most up to date displacement vector in load stepping procedure
* `displPrev`           displacement vector at end of previous load step
* `staticAnalysisSuccessful` boolean flag, true if load step was completed successfully
* `staticAnalysisComplete`   boolean flag, true if analysis is complete
"""
function updateLoadStep(iterationCount,loadStepParams,loadStep,loadStepPrev,loadStepCount,displPrev,displ)

    if loadStepParams.adaptiveLoadSteppingFlag #for adaptive load stepping option
        #check if maximum number of load steps has been exceeded.
        if loadStepCount > loadStepParams.maxNumLoadSteps
            error("Maximum number of load steps exceeded. Exiting.")
        end
        #calculate new loadstep adaptively
        loadStep,loadStepPrev,staticAnalysisSuccessful,staticAnalysisComplete = adaptiveLoadStepping(iterationCount,loadStepParams,loadStep,loadStepPrev)

    else #for prescribed load stepping option
        if iterationCount<loadStepParams.maxIterations #see if previous load step was successful
            staticAnalysisSuccessful= true
            if loadStep == 1.0 #if load step = 1.0, analysis is complete
                staticAnalysisComplete = true
            else
                staticAnalysisComplete = false
                loadStepPrev = loadStep
                loadStep = loadStepParams.prescribedLoadStep[loadStepCount]
                println("Prescribed load step: $loadStep")
            end
        else
            staticAnalysisSuccessful = false
        end

    end
    if staticAnalysisSuccessful
        displPrev = displ #update displacementPrev variable if previous load step was successful.
    else
        displ = displPrev #reset to displ vector to that at start of load step if load step was unsuccessful

        if !loadStepParams.adaptiveLoadSteppingFlag #if
            error("Maximum number of iterations exceeded in prescribed loadstep profile. Exiting.")
        end
    end
    return loadStep,loadStepPrev,displ,displPrev,staticAnalysisSuccessful,staticAnalysisComplete
end

"""
Internal, performs updates a loadstep adaptively, see ?updateLoadStep
"""
function adaptiveLoadStepping(iterationCount,loadStepParams,loadStep,loadStepPrev)

    staticAnalysisComplete = false
    loadStepOrig = loadStep

    if iterationCount>=loadStepParams.maxIterations                   #check for exceeding max iterations
        msgId = 1 #corresponds to a message  saying load step was unsuccessful and is being reduced
        staticAnalysisSuccessful = false
    else
        if loadStep == 1.0
            msgId = 2 #corresponds to a message saying loadstep was successful and load stepping is finished
        else
            msgId = 3 #corresponds to a message saying loads tep was successful and analysis is proceeding to next load step
        end
        staticAnalysisSuccessful=true
        if loadStep == 1.0
            staticAnalysisComplete = true
        end
        loadStepPrev = loadStep #update previous load step and end loadstep
        loadStep = 1.0
    end

    loadStepCopy = loadStep #make copy of load step for later checks
    loadStep = 0.5*(loadStep + loadStepPrev) #update load step


    if abs(loadStep-loadStepPrev)< loadStepParams.minLoadStepDelta #enforces delta load step is not below the minimum specified value
        if loadStep<loadStepPrev
            loadStep = loadStep - loadStepParams.minLoadStepDelta
        else
            loadStep = loadStep + loadStepParams.minLoadStepDelta
        end
    end

    if loadStep < loadStepParams.minLoadStep  #check that load step is not below the minimum specified load step
        loadStep = loadStepParams.minLoadStep
        if loadStepCopy == loadStep
            error("Minimum load step reached. Exiting.")
        end
    elseif loadStep > 1.0 #if load step has extended beyond 1.0, set to 1.0
        loadStep = 1.0
    end

    #print loadstepping message to command line
    if msgId == 1
        println("Max iterations exceeded for loadstep ($loadStepOrig) Reducing load step to $loadStep")
    elseif msgId == 2
        println("Nonlinear iteration successful for loadstep ($loadStep)  Nonlinear static analysis complete.")
    elseif msgId == 3
        println("Nonlinear iteration successful for loadstep ($loadStepOrig)  Increasing loadstep size to $loadStep")
    end

    return loadStep,loadStepPrev,staticAnalysisSuccessful,staticAnalysisComplete
end
