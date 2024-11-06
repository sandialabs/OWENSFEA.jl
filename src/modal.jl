"""
frequencies = autoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    minRPM = 0.0,
    maxRPM = 40.0,
    NRPM = 9, # int
    VTKsavename = nothing,
    saveModes = [1,3,5], #must be int
    saveRPM = [1,3,5], #must be int
    mode_scaling = 500.0,
    )

Automated Campbell Diagram Generator, this function runs the model with centrifugal stiffening for the specified RPM levels.  
If FEAinputs.analysisType == "GX" and VTKsavename are specified, it will output paraview mode shape files at the specified save name. 

#Inputs
* `FEAinputs::OWENSFEA.FEAModel`: The FEA modeling options
* `mymesh::OWENSFEA.Mesh`: a previously generated turbine mesh
* `myel::OWENSFEA.El`: the element properties associated with that mesh
* `system::GXBeam.System`: the converted GXBeam system from the mesh and el
* `assembly::GXBeam.AssemblyState`: the converted GXBeam assembly from the mesh and el
* `sections::Array{Float64, 3}`: the 3D point cloud to be converted to VTK format
* `minRPM::Float64`: minimum RPM to be run, e.x. 0.0
* `maxRPM::Float64`: maximum RPM to be run e.x. 40.0
* `NRPM::Int64`: number of linear discretizations of RPM e.x. 9 must be int
* `VTKsavename::string`: filename (with path if desired) of the VTK outputs if GX.  Set to "nothing" to not save.
* `saveModes::Array{Int64}`: The modes to save in the VTK outputs e.x. [1,3,5] must be int
* `saveRPM::Array{Int64}`: The RPMs to save in the VTK outputs e.x. [1,3,5] must be int
* `mode_scaling::Float64`: The mode scaling in the VTK outputs e.x. 500.0

#Outputs
* `frequency::Array{Float64}`: The output modal frequencies
"""
function autoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    minRPM = 0.0,
    maxRPM = 40.0,
    NRPM = 9, # int
    VTKsavename = nothing,
    saveModes = [1,3,5], #must be int
    saveRPM = [1,3,5], #must be int
    mode_scaling = 500.0,
    )

    rotSpdArrayRPM = LinRange(minRPM,maxRPM,NRPM)

    # Campbell Diagram generation
    if FEAinputs.analysisType!="GX"
        rotorSpeedArrayHZ = rotSpdArrayRPM ./ 60.0
        freq = zeros(length(rotorSpeedArrayHZ),FEAinputs.numModes)
        for irpm=1:length(rotorSpeedArrayHZ)
            println("$irpm of $(length(rotorSpeedArrayHZ))")
            rotorSpeed = rotorSpeedArrayHZ[irpm]
            local Omega = rotorSpeed
            displ = zeros(mymesh.numNodes*6)
            OmegaStart = Omega
            
            freqtemp,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,
            U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,displ=modal(FEAinputs,
            mymesh,myel;displ,Omega,OmegaStart,returnDynMatrices=false)
            
            freq[irpm,:] = freqtemp[1:FEAinputs.numModes]
        end
        
        return freq
    else
        ########################################
        ############ GXBeam ####################
        ########################################

        prescribed_conditions  = setPrescribedConditions(mymesh;pBC=FEAinputs.BC.pBC)#,Fexternal,ForceDof)
        
        # prescribed_conditions = Dict(
        #     # fixed base
        #     1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        #     # fixed top, but free to rotate around z-axis
        #     top_idx => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0),
        # )

        # --- Perform Analysis --- #

        # gravity vector
        gravity = [0, 0, -9.81] #TODO: from FEAinputs

        # number of modes
        nmode = FEAinputs.numModes

        # number of eigenvalues
        nev = 2*nmode

        # storage for results
        freq2 = zeros(length(rotSpdArrayRPM), nmode)
        λ_save = Array{Array{ComplexF64,1},1}(undef,length(rotSpdArrayRPM))
        eigenstates_save = Array{Array{GXBeam.AssemblyState,1},1}(undef,length(rotSpdArrayRPM))
        Up = []
        # perform an analysis for each rotation rate
        for (irpm,rpm) in enumerate(rotSpdArrayRPM)
            println("$irpm of $(length(rotSpdArrayRPM))")
            # set turbine rotation
            angular_velocity = [0, 0, rpm*(2*pi)/60]

            # eigenvalues and (right) eigenvectors
            system, λ, V, converged = GXBeam.eigenvalue_analysis!(system, assembly;
                prescribed_conditions = prescribed_conditions,
                angular_velocity = angular_velocity,
                gravity = gravity,
                nev = nev
                )

            # check convergence
            @assert converged

            if irpm > 1
                # construct correlation matrix
                C = Up*system.M*V

                # correlate eigenmodes
                perm, corruption = GXBeam.correlate_eigenmodes(C)

                # re-arrange eigenvalues
                λ = λ[perm]

                # update left eigenvector matrix
                Up = GXBeam.left_eigenvectors(system, λ, V)
                Up = Up[perm,:]
            else
                # update left eigenvector matrix
                Up = GXBeam.left_eigenvectors(system, λ, V)
            end

            # save frequencies
            freq2[irpm,:] = [imag(λ[k])/(2*pi) for k = 1:2:nev]

            state = GXBeam.AssemblyState(system, assembly;
            prescribed_conditions = prescribed_conditions)
            eigenstates = [GXBeam.AssemblyState(V[:,k],system, assembly;
            prescribed_conditions = prescribed_conditions) for k = 1:nev]
            λ_save[irpm] = λ
            eigenstates_save[irpm] = eigenstates
        end

        if !isnothing(VTKsavename) #TODO: map the OWENS state into the gx state for filesaving
            state = GXBeam.AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)
    
            try #this should error if someone on windows uses backslash '\'
                lastforwardslash = findlast(x->x=='/',VTKsavename)
                filepath = VTKsavename[1:lastforwardslash-1]
                if !isdir(filepath)
                    mkdir(filepath)
                end
            catch
                @info "Please manually create the directory to house $VTKsavename"
            end
            for isaveRPM in saveRPM
                for isavemode in saveModes
                    GXBeam.write_vtk("$(VTKsavename)_RPM$(rotSpdArrayRPM[isaveRPM])_Mode$(isavemode)", assembly, state, 
                        λ_save[isaveRPM][isavemode], eigenstates_save[isaveRPM][isavemode]; sections,mode_scaling)
                end
            end
        end

        return freq2
    end
end

"""
    modal(feamodel,mesh,el;Omega=0.0,displ=zeros(mesh.numNodes*6),OmegaStart=0.0,returnDynMatrices=false)

Modal analysis

# Inputs
* `feamodel::FEAModel`: see ?FEAModel
* `mesh::Mesh`: see ?Mesh
* `el::El`: see ?El
* `Omega::float`: Rotational rate in Hz
* `displ::Array{<:float}`: zeros(mesh.numNodes*6) initial (warm start) displacements for each dof
* `OmegaStart::float`: rotor speed (Hz) from previous analysis if stepping through various rotor speeds, may be useful in load stepping
* `returnDynMatrices::Bool`: Flag to save linearized K/C/M matrices for the design

# Outputs:
* `freq::Array{<:float}`: sorted modal frequencies (Hz)
* `damp::Array{<:float}`: sorted modal damping
* `imagCompSign::Array{<:float}`: sign of imaginary component of eigenvalues
* `U_x_0::Array{<:float}`: NnodesxNmodes in-phase mode shape x
* `U_y_0::Array{<:float}`: NnodesxNmodes in-phase mode shape y
* `U_z_0::Array{<:float}`: NnodesxNmodes in-phase mode shape z
* `theta_x_0::Array{<:float}`: NnodesxNmodes in-phase mode shape rotation about x
* `theta_y_0::Array{<:float}`: NnodesxNmodes in-phase mode shape rotation about y
* `theta_z_0::Array{<:float}`: NnodesxNmodes in-phase mode shape rotation about z
* `U_x_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape x
* `U_y_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape y
* `U_z_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape z
* `theta_x_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape rotation about x
* `theta_y_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape rotation about y
* `theta_z_90::Array{<:float}`: NnodesxNmodes out-of-phase mode shape rotation about z

"""
function modal(feamodel,mesh,el;Omega=0.0,displ=zeros(mesh.numNodes*6),OmegaStart=0.0,returnDynMatrices=false,elStorage=nothing,predef=nothing)

    if elStorage===nothing
        elStorage = initialElementCalculations(feamodel,el,mesh) #performs initial element calculations
    end
    # [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage) #calculate mass properties of structure

    #Do nonlinear iteration if needed
    if feamodel.spinUpOn
        feamodel.aeroElasticOn = false
        displ,_,staticAnalysisSuccessful=staticAnalysis(feamodel,mesh,el,displ,Omega,
        OmegaStart,elStorage) #performs static analysis about specified operating condition
    else
        staticAnalysisSuccessful = true
    end

    if staticAnalysisSuccessful
        freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,
        U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,eigVal,eigVec= linearAnalysisModal(feamodel,
        mesh,el,displ,Omega,elStorage;returnDynMatrices,predef) #performs modal analysis
    else
        error("Static analysis unsuccessful. Exiting")
    end
    return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,
    U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,displ,eigVal,eigVec
end

"""
Internal, see ?modal
"""
function  linearAnalysisModal(feamodel,mesh,el,displ,Omega,elStorage;returnDynMatrices=false,predef=nothing)

    feamodel.analysisType = "M" #Force type to align with the modal call
    elementOrder = feamodel.elementOrder  #extract element order from feamodel
    numNodesPerEl = elementOrder + 1  #do initialization
    numDOFPerNode = 6
    totalNumDOF = mesh.numNodes * numDOFPerNode
    countedNodes = []

    Kg = zeros(totalNumDOF,totalNumDOF)
    Mg = zeros(totalNumDOF,totalNumDOF)
    Cg = zeros(totalNumDOF,totalNumDOF)
    eldisp = zeros(numNodesPerEl*numDOFPerNode)

    timeInt = TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,
    displ,Omega,elStorage;Kg,Mg,Cg,countedNodes,predef)

    # #apply general 6x6  mass, damping, and stiffness matrices to nodes
    # Kg_all,Mg_all,Cg_all = applyGeneralConcentratedTerms(Kg,Mg,Cg,feamodel.nodalTerms.concStiffGen,feamodel.nodalTerms.concMassGen,feamodel.nodalTerms.concDampGen)

    #----------------------------------------------------------------------
    #APPLY CONSTRAINT
    Kg_con = applyConstraints(Kg,feamodel.jointTransform)  #modify global matrices for joint constraints using joint transform
    Mg_con = applyConstraints(Mg,feamodel.jointTransform)
    Cg_con = applyConstraints(Cg,feamodel.jointTransform)

    #APPLY BOUNDARY CONDITIONS
    KgTotal = applyBCModal(Kg_con,feamodel.BC.numpBC,feamodel.BC.map)     #apply boundary conditions to global matrices
    MgTotal = applyBCModal(Mg_con,feamodel.BC.numpBC,feamodel.BC.map)
    CgTotal = applyBCModal(Cg_con,feamodel.BC.numpBC,feamodel.BC.map)

    if Omega==0.0 #set eigensolver flag
        solveFlag = 2
    else
        solveFlag = 1
    end
    # eigVec,eigVal = eigSolve(MgTotal,CgTotal,KgTotal)#,... #eigensolve of global system
    if returnDynMatrices==true
        #save them to a file
        KgTotalU,_ = applyBC(Kg_con,zeros(length(Kg[:,1])),feamodel.BC,numDOFPerNode)
        MgTotalU,_ = applyBC(Mg_con,zeros(length(Mg[:,1])),feamodel.BC,numDOFPerNode)
        CgTotalU,_ = applyBC(Cg_con,zeros(length(Cg[:,1])),feamodel.BC,numDOFPerNode)

        filename = "./linearized_matrices.mat"
        println("Saving linearized matrices to: $filename")
        file = MAT.matopen(filename,"w")
        MAT.write(file,"Kg_all",Kg)
        MAT.write(file,"Mg_all",Mg)
        MAT.write(file,"Cg_all",Cg)
        MAT.write(file,"KgTotalU",KgTotalU)
        MAT.write(file,"MgTotalU",MgTotalU)
        MAT.write(file,"CgTotalU",CgTotalU)
        MAT.write(file,"KgTotalM",KgTotal)
        MAT.write(file,"MgTotalM",MgTotal)
        MAT.write(file,"CgTotalM",CgTotal)
        MAT.write(file,"jointTransform",feamodel.jointTransform)
        MAT.write(file,"BC",feamodel.BC)
        MAT.write(file,"mesh",mesh)
        close(file)
    end
    #constructs state space form (with mass matrix inverted)
    matwidth = length(MgTotal[:,1])
    sysMat = [zeros(matwidth,matwidth) 1.0*LinearAlgebra.I(matwidth)
    -MgTotal\KgTotal -MgTotal\CgTotal]
    nev = size(sysMat)[1]#min(size(sysMat)[1],feamodel.numModes)
    # eigVal, eigVec = ArnoldiMethod.partialeigen(ArnoldiMethod.partialschur(sysMat))# nev=min(nx,nev), which=ArnoldiMethod.LM())[1])
    eigVal, eigVec = ArnoldiMethod.partialeigen(ArnoldiMethod.partialschur(sysMat; nev, which=ArnoldiMethod.LM())[1])
    # println(maximum(abs.(sysMat*eigVec .- eigVec*(eigVal.*LinearAlgebra.I(length(eigVal))))))
    #feamodel.numModesToExtract,solveFlag)

    perm = sortperm(eigVal, by=(eigVal)->(abs(eigVal),imag(eigVal)), rev=false)
    eigVal .= eigVal[perm]
    eigVec .= eigVec[:,perm]

    #extract frequency, damping, mode shapes from eigenvalues and vectors
    len3 = length(eigVal)
    freq = zeros(len3)
    damp = zeros(len3)
    len1 = Int(length(displ)/numDOFPerNode)
    phase1 = zeros(len1,numDOFPerNode,len3)
    phase2 = zeros(len1,numDOFPerNode,len3)
    sortedModes0 = zeros(len1,numDOFPerNode,len3)
    sortedModes = complex.(sortedModes0,zero(sortedModes0))
    imagCompSign = zeros(len3)
    for i=1:len3
        freq[i],damp[i],phase1[:,:,i],phase2[:,:,i],sortedModes[:,:,i] = extractFreqDamp(eigVal[i],eigVec[:,i],numDOFPerNode,feamodel.jointTransform,feamodel.reducedDOFList,feamodel.BC,feamodel.analysisType)
        imagCompSign[i] = sign(imag(eigVal[i]))
    end

    # #write output
    if feamodel.analysisType !="FA"
        freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90 = ModalOutput(freq,damp,phase1,phase2,imagCompSign,feamodel.outFilename)
    end

    return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,eigVal,eigVec

end

"""
    applyBCModal(K,BC,numDofPerNode)

Internal, applies boundary conditions to a system matrix for modal analysis

# Inputs
* `K`:             assembled global system matrix
* `BC`:            struct of boundary condition information
* `numDofPerNode`: number of degrees of freedom per node

# Outputs:
* `Knew` global system matrix with boundary conditions
"""
function applyBCModal(K,numpBC,bcMap)

    numEq=length(bcMap)
    # indVec = zeros(numEq-numpBC,1) #can't pre-initialize...puts zeros in map
    # that causes problems when creating Knew

    index = 1
    indVec = zeros(Int,sum(bcMap .!= -1))
    for i=1:numEq
        if bcMap[i] != -1
            indVec[index] = i
            index = index +1
        end
    end

    if numpBC > 0
        Knew = K[indVec,indVec]
    else
        Knew = K
    end
    return Knew
end

"""
    assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)

Internal, assembles the element matrix into the global system of equations

# Inputs
* `Ke`:             element matrix
* `conn`:           element connectivity
* `numNodesPerEl`:  number of nodes per element
* `numDofPerNode`:  number of degrees of freedom per node
* `Kg`:             global system matrix

# Outputs:
* `Kg`:             global system matrix with assembled element
"""
function assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)

    count = 1
    dofList = zeros(Int,numNodesPerEl*numDOFPerNode)
    for i=1:numNodesPerEl
        for j=1:numDOFPerNode
            dofList[count] = (conn[i]-1)*numDOFPerNode + j
            count = count +1
        end
    end

    Kg[dofList,dofList] = Kg[dofList,dofList] + Ke
    # numDOFPerEl = length(dofList)
    # #Assemble element i into global system
    #         for j=1:numDOFPerEl
    #             J = dofList(j)
    #             for m=1:numDOFPerEl
    #                 M = dofList(m)
    #                 Kg(J,M) = Kg(J,M) + Ke(j,m)
    #             end
    #         end

    return Kg

end

function assemblyMatrixOnly!(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)

    count = 1
    dofList = zeros(Int,numNodesPerEl*numDOFPerNode)
    for i=1:numNodesPerEl
        for j=1:numDOFPerNode
            dofList[count] = (conn[i]-1)*numDOFPerNode + j
            count = count +1
        end
    end

    Kg[dofList,dofList] = Kg[dofList,dofList] + Ke
end

"""
    extractFreqDamp(val,vec,numDOFPerNode,jointTransform,reducedDOFList,BC,analysisType)

Internal, calculates the eigenvalues and vectors of a structural dynamic system

# Inputs
* `val`:             eigenvalue
* `vec`:             eigenvector
* `numDOFPerNode`:   number of degrees of freedom per node
* `jointTransform`:  joint transformation matrix from reduced to full DOF list
* `reducedDOFList`:  listing of reduced DOFs
* `BC`:              boundary condition object containing boundary condition info
* `analysisType`:    analysis type

# Outputs:
* `freq`:        modal frequency
* `damp`:        modal damping
* `phase1`:      in phase mode shape (real part of mode shape)
* `phase2`:      out of phase mode shape (imaginary part of mode shape)
* `sortedModes`: total, complex mode shape
"""
function extractFreqDamp(val,vec,numDOFPerNode,jointTransform,reducedDOFList,BC,analysisType)

    freq = abs(imag(val))/(2*pi)   #damped frequency for state space system
    damp = -real(val)/abs(imag(val)) #modal damping

    if (abs(imag(val)) < 1.0e-4)     #if imaginary part of eigenvalue is numeric zero treat as spring-mass system
        freq = sqrt(abs(real(val)))/(2*pi)
        damp = 0.0
    end

    # if (~strcmp(analysisType,'FA'))  #for all but automated flutter analysis
    #          [len,numModeShapes] = size(vec)
    dispReduc = constructReducedDispVecFromEigVec(vec,reducedDOFList,BC) #construct mode shape vector with boundary conditinos
    dispOrig = jointTransform*dispReduc #transform from reduced DOF list to full DOF list
    lenOrig=length(dispOrig)

    sortedModes0=zeros(Int(lenOrig/numDOFPerNode),numDOFPerNode)
    sortedModes = complex.(sortedModes0,0)
    for i=1:Int(lenOrig/numDOFPerNode)     #construct matrix of nodal DOF values from full DOF eigenvector
        for j=1:numDOFPerNode
            sortedModes[i,j] = dispOrig[(i-1)*6+j]
        end
    end

    phase1 = real(sortedModes)  #phase 1 is real part of modeshape (0 deg in phase)
    phase2 = imag(sortedModes)  #phase 2 is imag part of modeshape (90 deg out of phase)

    max1=maximum(abs.(phase1)) #find maximum values for modeshape normalization
    max2=maximum(abs.(phase2))


    if (max1>max2)
        maxOverall = max1
    else
        maxOverall = max2
    end

    if (maxOverall == 0)
        maxOverall = 1
    end

    phase1 = phase1./maxOverall  #normalize modeshapes
    phase2 = phase2./maxOverall

    if (abs(minimum(phase1))+1)<1.0e-4 || (abs(minimum(phase2))+1)<1.0e-4
        phase1 = -1*phase1
        phase2 = -1*phase2

    end

    return freq,damp,phase1,phase2,sortedModes
end

"""
    constructReducedDispVecFromEigVec(vec1,reducedDOFList,BC)

Internal, This function takes the original mode shape and modifies it to account for boundary conditions

# Inputs
* `vec1`:
* `reducedDOFList`:
* `BC`:

# Outputs:
* `vec1Red`:
"""
function constructReducedDispVecFromEigVec(vec1,reducedDOFList,BC)
    #This function takes the original mode shape and modifies it to
    #account for boundary conditions
    bcdoflist=zeros(1,BC.numpBC)
    #form pBC DOF list
    for i=1:BC.numpBC
        bcnodenum = BC.pBC[i,1]
        bcdofnum = BC.pBC[i,2]
        bcdoflist[i] = (bcnodenum-1)*6 + bcdofnum
    end

    index = 1
    len = length(vec1)
    vec1Red0 = zeros(length(reducedDOFList))
    vec1Red = complex.(vec1Red0,0)
    for i=1:length(reducedDOFList)
        if isnothing(findfirst(x->x==reducedDOFList[i],bcdoflist))#(!ismember(reducedDOFList[i],bcdoflist))
            vec1Red[i] = vec1[Int(len/2+index),1]
            index = index + 1
        else
            vec1Red[i] = 0
        end
    end
    return vec1Red
end
