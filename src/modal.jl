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
function modal(feamodel,mesh,el;Omega=0.0,displ=zeros(mesh.numNodes*6),OmegaStart=0.0,returnDynMatrices=false)

    elStorage = initialElementCalculations(feamodel,el,mesh) #performs initial element calculations

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
        U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90= linearAnalysisModal(feamodel,
        mesh,el,displ,Omega,elStorage;returnDynMatrices) #performs modal analysis
    else
        error("Static analysis unsuccessful. Exiting")
    end
    return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,
    U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,displ
end

"""
Internal, see ?modal
"""
function  linearAnalysisModal(feamodel,mesh,el,displ,Omega,elStorage;returnDynMatrices=false)

    feamodel.analysisType = "M" #Force type to align with the modal call
    elementOrder = feamodel.elementOrder  #extract element order from feamodel
    numNodesPerEl = elementOrder + 1  #do initialization
    numDOFPerNode = 6
    totalNumDOF = mesh.numNodes * numDOFPerNode

    Kg = zeros(totalNumDOF,totalNumDOF)
    Mg = zeros(totalNumDOF,totalNumDOF)
    Cg = zeros(totalNumDOF,totalNumDOF)
    eldisp = zeros(numNodesPerEl*numDOFPerNode)

    nodalTerms,timeInt = TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,
    displ,Omega,elStorage;Kg,Mg,Cg)

    #apply general 6x6  mass, damping, and stiffness matrices to nodes
    Kg_all,Mg_all,Cg_all = applyGeneralConcentratedTerms(Kg,Mg,Cg,feamodel.nodalTerms.concStiffGen,feamodel.nodalTerms.concMassGen,feamodel.nodalTerms.concDampGen)

    #----------------------------------------------------------------------
    #APPLY CONSTRAINT
    Kg = applyConstraints(Kg_all,feamodel.jointTransform)  #modify global matrices for joint constraints using joint transform
    Mg = applyConstraints(Mg_all,feamodel.jointTransform)
    Cg = applyConstraints(Cg_all,feamodel.jointTransform)

    #APPLY BOUNDARY CONDITIONS
    KgTotal = applyBCModal(Kg,feamodel.BC.numpBC,feamodel.BC.map)     #apply boundary conditions to global matrices
    MgTotal = applyBCModal(Mg,feamodel.BC.numpBC,feamodel.BC.map)
    CgTotal = applyBCModal(Cg,feamodel.BC.numpBC,feamodel.BC.map)

    if Omega==0.0 #set eigensolver flag
        solveFlag = 2
    else
        solveFlag = 1
    end
    # eigVec,eigVal = eigSolve(MgTotal,CgTotal,KgTotal)#,... #eigensolve of global system
    if returnDynMatrices==true
        #save them to a file
        KgTotalU,_ = applyBC(Kg,zeros(length(Kg[:,1])),feamodel.BC,numDOFPerNode)
        MgTotalU,_ = applyBC(Mg,zeros(length(Mg[:,1])),feamodel.BC,numDOFPerNode)
        CgTotalU,_ = applyBC(Cg,zeros(length(Cg[:,1])),feamodel.BC,numDOFPerNode)

        filename = "./linearized_matrices.mat"
        println("Saving linearized matrices to: $filename")
        file = MAT.matopen(filename,"w")
        MAT.write(file,"Kg_all",Kg_all)
        MAT.write(file,"Mg_all",Mg_all)
        MAT.write(file,"Cg_all",Cg_all)
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

    return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90

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
