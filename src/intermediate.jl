
function TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,dispData,Omega,elStorage;
    Kg=zeros(1),Mg=zero(Kg),Cg=zero(Kg),Fg=zero(Kg[:,1]),
    eldisp_sm1=zero(eldisp),eldispdot=zero(eldisp),eldispddot=zero(eldisp),
    eldispiter=zero(eldisp),rbData=zeros(9),CN2H=1.0*LinearAlgebra.I(3),delta_t=0.0,
    OmegaDot=0.0,displ_im1=zero(eldisp),displdot_im1=zero(eldisp),timeInt=nothing,
    displddot_im1=zero(eldisp),iterationCount=0,postprocess=false,elementNumber=1,
    dispOld=nothing,loadStep = 1.0,loadStepPrev = 1.0,predef=nothing,countedNodes = [])

    x = mesh.x
    y = mesh.y
    z = mesh.z
    conn = Int.(mesh.conn)

    elementOrder = feamodel.elementOrder
    numNodesPerEl = elementOrder + 1
    numDOFPerNode = 6
    totalNumDOF = mesh.numNodes * numDOFPerNode
    analysisType = feamodel.analysisType
    modalFlag = true
    aeroForceOn = false

    nlParams = feamodel.nlParams
    gravityOn = feamodel.gravityOn
    RayleighAlpha = feamodel.RayleighAlpha
    RayleighBeta = feamodel.RayleighBeta
    accelVec = rbData[1:3] #TODO: check that if analysis type is not unsteady, these three vecs should be zero
    omegaVec = rbData[4:6]
    omegaDotVec = rbData[7:9]
    airDensity = feamodel.airDensity
    maxNumLoadSteps = nlParams.maxNumLoadSteps
    MAXIT = nlParams.maxIterations
    tolerance = nlParams.tolerance

    #------- intitialization -----------------
    if analysisType=="TNB"
        if timeInt == nothing
            #------ newmark integration parameters ---------
            alpha = 0.5
            gamma = 0.5
            beta = 0.5*gamma

            a1 = alpha*delta_t
            a2 = (1.0-alpha)*delta_t
            a3 = 1.0/(beta*delta_t*delta_t)
            a3 = a3
            a4 = a3*delta_t
            a5 = 1.0/gamma-1.0
            a6 = alpha/(beta*delta_t)
            a7 = alpha/beta - 1.0
            a8 = delta_t*(alpha/gamma-1.0)

            timeInt = TimeInt(delta_t,a1,a2,a3,a4,a5,a6,a7,a8)
        end

        disp_s = dispData.displ_s
        dispdot_s = dispData.displdot_s
        dispddot_s = dispData.displddot_s

        aeroElasticOn = false
        iterationType = feamodel.nlParams.iterationType
        useDisp = feamodel.nlOn
        preStress = false

    elseif analysisType=="TD"
        if timeInt == nothing
            #------ dean integration parameters -------------
            alpha = 0.25

            a1 = alpha*delta_t^2
            a2 = (1-2*alpha)*delta_t^2
            a3 = delta_t/2.0
            a4 = delta_t*delta_t

            timeInt = TimeInt(delta_t,a1,a2,a3,a4,0.0,0.0,0.0,0.0)
        end

        disp_s = dispData.displ_s
        disp_sm1 = dispData.displ_sm1

        aeroElasticOn = false
        iterationType = feamodel.nlParams.iterationType
        useDisp = feamodel.nlOn
        preStress = false
        #-------------------------------------------
    elseif analysisType=="M"
        if timeInt == nothing
            timeInt = TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        end
        aeroElasticOn = feamodel.aeroElasticOn
        iterationType = "NONE"
        useDisp = false
        preStress = true
        displ = dispData
    elseif analysisType=="S"
        if timeInt == nothing
            timeInt = TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        end
        aeroElasticOn = false
        if feamodel.nlOn
            iterationType = nlParams.iterationType  #define nonlinear iteration type
        else
            iterationType = "LINEAR"
        end
        useDisp = feamodel.nlOn
        preStress = false
        displ = dispData
        analysisType = "M" #steady and modal are equivalent for timoskenko?
        modalFlag = false
    elseif analysisType=="ROM"
        if timeInt == nothing
            timeInt = TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        end
        useDisp = true
        preStress = false
        iterationType = "LINEAR" #linear analysis option
        aeroElasticOn = feamodel.aeroElasticOn
        if postprocess
            analysisType = "TNB"
            disp_s = dispData.displ_s
            dispdot_s = dispData.displdot_s
            dispddot_s = dispData.displddot_s
        else
            displ = dispData
            airDensity = 0.0
        end
    else
        error("analysis type not supported, choose another")
    end

    if postprocess
        numEl = 1
        iterationType = "DI"
        useDisp = feamodel.nlOn
        preStress = false
    else
        numEl = mesh.numEl
    end

    elx=zeros(numNodesPerEl)
    ely=zeros(numNodesPerEl)
    elz=zeros(numNodesPerEl)

    #-------------------------------------------

    #---- element  calculation and assembly ----------------------------------
    for i=1:numEl
        if postprocess
            i = elementNumber
        end
        #Calculate Ke and Fe for element i
        index = 1                           #initialize element data
        xloc = [0.0 el.elLen[i]]
        sectionProps = el.props[i]
        sweepAngle = el.psi[i]
        coneAngle = el.theta[i]
        rollAngle = el.roll[i]
        aeroSweepAngle = 0.0
        if ((iterationCount == 0 || analysisType=="M") && postprocess == false) && analysisType!="ROM"
            firstIteration = true
        else
            firstIteration = false
        end

        if feamodel.analysisType == "S"
            if iterationCount >=1
                gamm = 0.0     #option for acceleration of iterative procedure (gamma = 0 or gamma=0.5 are typical)
                disp = dispOld.*gamm + displ.*(1-gamm)
                firstIteration = false
            else
                firstIteration = true
            end
        end

        for j=1:numNodesPerEl
            #get element cooridnates
            elx[j] = x[conn[i,j]]
            ely[j] = y[conn[i,j]]
            elz[j] = z[conn[i,j]]
            #get element nodal displacements at s and s-1 time step
            for k=1:numDOFPerNode
                if analysisType=="TD"
                    eldisp[index] = disp_s[(conn[i,j]-1)*numDOFPerNode + k]
                    eldisp_sm1[index] = disp_sm1[(conn[i,j]-1)*numDOFPerNode + k]
                    eldispiter[index] = displ_iter[(conn[i,j]-1)*numDOFPerNode + k]
                elseif analysisType=="TNB"
                    eldispiter[index] = displ_im1[(conn[i,j]-1)*numDOFPerNode + k]
                    if (iterationType=="NR")
                        eldisp[index] = displ_im1[(conn[i,j]-1)*numDOFPerNode + k]
                        eldispdot[index] = displdot_im1[(conn[i,j]-1)*numDOFPerNode + k]
                        eldispddot[index] = displddot_im1[(conn[i,j]-1)*numDOFPerNode + k]
                    elseif (iterationType=="DI"||iterationType=="LINEAR")
                        eldisp[index] = disp_s[(conn[i,j]-1)*numDOFPerNode + k]
                        eldispdot[index] = dispdot_s[(conn[i,j]-1)*numDOFPerNode + k]
                        eldispddot[index] = dispddot_s[(conn[i,j]-1)*numDOFPerNode + k]
                    end
                elseif analysisType=="M" || analysisType=="ROM"
                    eldisp[index] = displ[(conn[i,j]-1)*numDOFPerNode + k]
                end
                index = index + 1
            end
        end

        #get concentrated terms associated with elemet
        concLoad
        concStiff
        concMass
        concDamp
        # concStiff, concMass, concDamp, concLoad = applyConcentratedTerms(feamodel.nodalTerms.concStiff, feamodel.nodalTerms.concMass, feamodel.nodalTerms.concDamp, feamodel.nodalTerms.concLoad, feamodel.joint, numNodes, numDOFPerNode)
        concStiff, concMass, concDamp, concLoad, countedNodes = getElementConcTerms!(feamodel.nodalTerms.concStiff, feamodel.nodalTerms.concMass, feamodel.nodalTerms.concDamp, feamodel.nodalTerms.concLoad, conn[i,:], numDOFPerNode, countedNodes)
        if el.rotationalEffects[i]!=1
            Omega = 0.0
            OmegaDot = 0.0
    		omegaVec = zeros(3)
    		omegaDotVec = zeros(3)
        end

        if aeroElasticOn #TODO: why are we doing this if it is done again within timoshenko?
            freq = feamodel.guessFreq*2.0*pi     #set guess frequency if aeroelastic analysis
        else
            freq = 0.0 #Declare variable on all execution paths
        end

        if analysisType == "TD" || analysisType == "TNB" || postprocess
            freq = 0.0
        end

        if !(analysisType=="TNB" || analysisType=="ROM") #TODO: TD?
            eldispiter = eldisp
        end

        elInput = ElInput(elementOrder,modalFlag,timeInt,xloc,sectionProps,sweepAngle,
        coneAngle,rollAngle,aeroSweepAngle,iterationType,useDisp,preStress,aeroElasticOn,
        aeroForceOn,loadStepPrev,loadStep,maxNumLoadSteps,MAXIT,tolerance,analysisType,
        eldisp,eldispdot,eldispddot,eldispiter,concMass,concStiff,concDamp,concLoad,eldisp_sm1,elx,ely,elz,
        gravityOn,RayleighAlpha,RayleighBeta,accelVec,omegaVec,omegaDotVec,Omega,
        OmegaDot,CN2H,airDensity,freq,firstIteration)

        elOutput = calculateTimoshenkoElementNL(elInput,elStorage[i];predef) #calculate timoshenko element
        if postprocess
                 #post process for reaction force
                FhatEl1PP = elOutput.Ke*eldispiter
                if analysisType=="TD"
                    denom = timeInt.a4
                else
                    denom = 1.0
                end

                return (FhatEl1PP - elOutput.FhatLessConc)./denom

        elseif analysisType == "TD" || analysisType == "TNB" || feamodel.analysisType == "S"
            assembly!(elOutput.Ke,elOutput.Fe,conn[i,:],numNodesPerEl,numDOFPerNode,Kg,Fg) #assemble element stiffness matrix and force vector
        elseif analysisType == "M"
            assemblyMatrixOnly!(elOutput.Ke,conn[i,:],numNodesPerEl,numDOFPerNode,Kg) #assemble element into global stiffness matrix
            assemblyMatrixOnly!(elOutput.Me,conn[i,:],numNodesPerEl,numDOFPerNode,Mg) #assemble element into global mass matrix
            assemblyMatrixOnly!(elOutput.Ce,conn[i,:],numNodesPerEl,numDOFPerNode,Cg) #assemble element into global damping matrix
        elseif analysisType == "ROM"
            assembly(elOutput.Ke,elOutput.Fe,conn[i,:],numNodesPerEl,numDOFPerNode,Kg,Fg) #assemble element Kg and Fg
            assemblyMatrixOnly(elOutput.Ce,conn[i,:],numNodesPerEl,numDOFPerNode,Cg) #assemble Mg
            assemblyMatrixOnly(elOutput.Me,conn[i,:],numNodesPerEl,numDOFPerNode,Mg) #assemble Mg
        end
    end #for
    return nodalTerms,timeInt
end
