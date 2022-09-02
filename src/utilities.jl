"""

   calculateStructureMassProps(elStorage)

This function caclulates structural mass properties of the finite element mesh (mass, moment of inertia, mass center) about the origin of the mesh coordinate system.

#Input
* `elStorage::ElStorage`    see ?ElStorage, object containing arrays of stored element information

#Output
* `structureMass::float`       mass of structure
* `structureMOI::float`        moment of inertia tensor of structgure
* `structureMassCenter::float` center of mass of structure
"""
function calculateStructureMassProps(elStorage)

    numElements = length(elStorage) #get number of elements

    structureMass = 0.0 #initialize structure mass and moment of inertia
    structureMOI = zeros(3,3)
    temp = zeros(3,1)
    for i=1:numElements #sum over elemetns contribution to mass and moment of inertia
        structureMass += elStorage[i].mel
        structureMOI .+= elStorage[i].moiel
        temp .+= elStorage[i].xmel
    end

    structureMassCenter = temp./structureMass #calculate mass center

    #modify moment of inertia to be about structure mass center
    x = structureMassCenter[1]
    y = structureMassCenter[2]
    z = structureMassCenter[3]

    structureMOI = structureMOI - structureMass*[(y^2+z^2) -x*y -x*z
    -x*y (x^2+z^2) -y*z
    -x*z -y*z (x^2+y^2)]


    return structureMass,structureMOI,structureMassCenter

end

"""

    calculateLambda(theta1,theta2,theta3)

This function calculates a transformation matrix to transform the element degree of freedom vector (12 DOFs) from the hub frame to the element frame. The transformation matrix is constructed via the direction cosine matrices of a 3-2-1 Euler rotation sequence.

#Input
*`theta1::float`:   angle (rad) of rotation for 1st rotation of 3-2-1 sequence
*`theta2::float`:   angle (rad) of rotation for 2nd rotation of 3-2-1 sequence
*`theta3::float`:   angle (rad) of rotation for 3rd rotation of 3-2-1 sequence

#Output
*`lambda::Array{<:float}`:   12 x 12 transformation matrix
"""
function calculateLambda(theta1,theta2,theta3)

    ct1 = cos(theta1); st1=sin(theta1);
    ct2 = cos(theta2); st2=sin(theta2);
    ct3 = cos(theta3); st3=sin(theta3);

    fac1 = st3*st2
    fac2 = ct3*st2
    dcm = [ct2*ct1           ct2*st1          -st2
    fac1*ct1-ct3*st1  fac1*st1+ct3*ct1  st3*ct2
    fac2*ct1+st3*st1  fac2*st1-st3*ct1  ct3*ct2]


    lambda = zeros(12,12)
    lambda[1:3,1:3] = dcm
    lambda[4:6,4:6] = dcm
    lambda[7:9,7:9] = dcm
    lambda[10:12,10:12] = dcm

    return lambda

end

"""
Internal, creates a constraint transformation matrix for a single joint. Tda is this matrix, dDOF contains a listing of dependent global DOFs associated with this joint, and aDOF contains a listing of active global DOFs associated with this joint.
"""
function createTda(jointType,slaveNodeNum,masterNodeNum,psi,theta,joint)

    if (jointType == 4 && (abs(abs(theta)-90)<1.0e-3 || (abs(abs(theta)-270)<1.0e-3) ))
        theta = 0.0
        jointType = 3
    end

    #calculate transformation matrix from hub frame to joint frame
    Lambda = calculateLambda(psi*pi/180.0,theta*pi/180.0,0.0)

    #Tda is a local mapping of dependent DOFs to active DOFs at a node
    # u_d = Tda * u_a
    # such that u_d is a list of local dependent slave DOFs at a jont and
    # u_a is a list of local dependent slave DOFs at a joint.


    if (jointType == 0) #for weld/fixed joint type
        activeDof0 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof0 = [1 2 3 4 5 6] #slave DOF list at joint

        #determine local active DOFs associated with slave node
        slaveActiveDof0 = determineActiveDofsFromSlaveNode(slaveDof0,6)

        Rda0 = -1.0*LinearAlgebra.I(6) #from constraint equation for fixed joint
        Rdd0 = 1.0*LinearAlgebra.I(6)

        Tda,dDOF,aDOF = getNodeMaps(Rdd0,Rda0,masterNodeNum,slaveNodeNum,slaveDof0,activeDof0,slaveActiveDof0)

    elseif (jointType == 1) #for pinned joint type
        activeDof1 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof1 = [1 2 3] #slave DOF list at joint

        #determine local active DOFs associated with slave node
        slaveActiveDof1 = determineActiveDofsFromSlaveNode(slaveDof1,6)

        Rda1 = -[1.0*LinearAlgebra.I(3);; zeros(3,6)] #from constraint equation for pinned joint
        Rdd1 = 1.0*LinearAlgebra.I(3)

        Tda,dDOF,aDOF = getNodeMaps(Rdd1,Rda1,masterNodeNum,slaveNodeNum,slaveDof1,activeDof1,slaveActiveDof1)

    elseif (jointType == 2)     #hinge axis along localy "2" frame of joint

        if ((abs(abs(psi)-90))<1.0e-3 || (abs(abs(psi)-270))<1.0e-3)
            activeDof2 = [1 2 3 4 5 6]
            slaveDof2  = [1 2 3 5 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof2 = determineActiveDofsFromSlaveNode(slaveDof2,6)

            globalConstraintEqMatrix2 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3) zeros(3,3)
            zeros(2,3) -[0 1 0;0 0 1] zeros(2,3) [0 1 0;0 0 1]]

        else
            activeDof2 = [1 2 3 4 5 6]
            slaveDof2 = [1 2 3 4 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof2 = determineActiveDofsFromSlaveNode(slaveDof2,6)

            localConstraintEqMatrix2 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 0 1] zeros(2,3) [1 0 0;0 0 1]]
            globalConstraintEqMatrix2 = localConstraintEqMatrix2*Lambda
        end
        #extract Rda from globalConstraintEqMatrix2
        Rda2 = zeros(length(globalConstraintEqMatrix2[:,1]),length(activeDof2)+length(slaveActiveDof2))
        index = 1
        for i=1:length(activeDof2)
            ind = activeDof2[i]
            Rda2[:,index] = globalConstraintEqMatrix2[:,ind]
            index += 1
        end

        for i=1:length(slaveActiveDof2)
            ind = Int(slaveActiveDof2[i]+6)
            Rda2[:,index] = globalConstraintEqMatrix2[:,ind]
            index += 1
        end

        #extract Rdd from globalConstraintEqMatrix2
        Rdd2 = zeros(length(globalConstraintEqMatrix2[:,1]),length(slaveDof2))
        for i=1:length(slaveDof2)
            ind = Int(slaveDof2[i]+6)
            Rdd2[:,i] = globalConstraintEqMatrix2[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd2,Rda2,masterNodeNum,slaveNodeNum,slaveDof2,activeDof2,slaveActiveDof2)

    elseif (jointType == 3)     #hinge axis along local "1" frame of joint

        if ((abs(abs(theta)-90))<1.0e-3 || (abs(abs(theta)-270))<1.0e-3)
            activeDof3 = [1 2 3 4 5 6]
            slaveDof3  = [1 2 3 4 5]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            globalConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 1 0] zeros(2,3) [1 0 0;0 1 0]]

        elseif ((abs(abs(psi)-90))<1.0e-3 || (abs(abs(psi)-270))<1.0e-3)

            activeDof3 = [1 2 3 4 5 6]
            slaveDof3  = [1 2 3 4 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            globalConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 0 1] zeros(2,3) [1 0 0;0 0 1]]

        else

            activeDof3 = [1 2 3 4 5 6]
            slaveDof3 = [1 2 3 5 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            localConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[0 1 0;0 0 1] zeros(2,3) [0 1 0;0 0 1]]
            globalConstraintEqMatrix3 = localConstraintEqMatrix3*Lambda

        end

        #extract Rda from globalConstraintEqMatrix3
        Rda3 = zeros(length(globalConstraintEqMatrix3[:,1]),length(activeDof3)+length(slaveActiveDof3))
        index = 1
        for i=1:length(activeDof3)
            ind = activeDof3[i]
            Rda3[:,index] = globalConstraintEqMatrix3[:,ind]
            index += 1
        end

        for i=1:length(slaveActiveDof3)
            ind = slaveActiveDof3[i]+6
            Rda3[:,index] = globalConstraintEqMatrix3[:,ind]
            index += 1
        end

        #extract Rdd from globalConstraintEqMatrix3
        Rdd3 = zeros(length(globalConstraintEqMatrix3[:,1]),length(slaveDof3))
        for i=1:length(slaveDof3)
            ind = slaveDof3[i]+6
            Rdd3[:,i] = globalConstraintEqMatrix3[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd3,Rda3,masterNodeNum,slaveNodeNum,slaveDof3,activeDof3,slaveActiveDof3)

    elseif (jointType == 4)     #hinge axis along local "3" frame of joint

        activeDof4 = [1 2 3 4 5 6]
        slaveDof4  = [1 2 3 4 5]

        #determine local active DOFs associated with slave node
        slaveActiveDof4 = determineActiveDofsFromSlaveNode(slaveDof4,6)

        localConstraintEqMatrix4 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
        zeros(2,3) -[1 0 0;0 1 0] zeros(2,3) [1 0 0;0 1 0]]
        globalConstraintEqMatrix4 = localConstraintEqMatrix4*Lambda

        #extract Rda from globalConstraintEqMatrix4
        Rda4 = zeros(length(globalConstraintEqMatrix4[:,1]),length(activeDof4)+length(slaveActiveDof4))
        index = 1
        for i=1:length(activeDof4)
            ind = activeDof4[i]
            Rda4[:,index] = globalConstraintEqMatrix4[:,ind]
            index += 1
        end

        for i=1:length(slaveActiveDof4)
            ind = slaveActiveDof4[i]+6
            Rda4[:,index] = globalConstraintEqMatrix4[:,ind]
            index += 1
        end

        #extract Rdd from globalConstraintEqMatrix4
        Rdd4 = zeros(length(globalConstraintEqMatrix4[:,1]),length(slaveDof4))
        for i=1:length(slaveDof4)
            ind = slaveDof4[i]+6
            Rdd4[:,i] = globalConstraintEqMatrix4[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd4,Rda4,masterNodeNum,slaveNodeNum,slaveDof4,activeDof4,slaveActiveDof4)


    elseif (jointType == 5) #rigid bar constraint type
        activeDof5 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof5 = [1 2 3 4 5 6] #slave DOF list at joint
        #determine local active DOFs associated with slave node
        slaveActiveDof5 = determineActiveDofsFromSlaveNode(slaveDof5,6)

        Rdd5 = 1.0*LinearAlgebra.I(6);  #need to define lx,ly,lz, from mesh level
        lx = joint[5]
        ly = joint[6]
        lz = joint[7]

        Rda5 = -1.0*LinearAlgebra.I(6)
        Rda5[1:3,4:6] = [0 -lz ly;lz 0 -lx;-ly lx 0]

        Tda,dDOF,aDOF = getNodeMaps(Rdd5,Rda5,masterNodeNum,slaveNodeNum,slaveDof5,activeDof5,slaveActiveDof5);

    else
        error("Correct jointType not specified, should be 1, 2, 3, 4, or 5")
    end
    return Tda,dDOF,aDOF
end

"""
Internal, gets node mapping
"""
function getNodeMaps(Rdd,Rda,masterNodeNum,slaveNodeNum,slaveDof,activeDof,slaveActiveDof)

    if (abs(LinearAlgebra.det(Rdd)) < 1.0e-3)
        error("Singular joint transformation matrix. Exiting")
    end

    Tda = -Rdd\Rda #calculate Tda

    numSlaveDOFs = length(slaveDof) #get number of joint DOFs for this joint
    numActiveDOFsFromMasterNode = length(activeDof) #get number of active DOFs for this joint

    dDOF = zeros(Int,numSlaveDOFs) #initialize arrays
    aMap = zeros(Int,numActiveDOFsFromMasterNode,1)

    for i=1:numSlaveDOFs
        #get global DOF numbers of slave DOFs for this joint
        dDOF[i,1] = (slaveNodeNum-1)*6 + slaveDof[i]
    end

    for i=1:numActiveDOFsFromMasterNode
        #get global DOF numbers of active DOFs for this joint from master nodes
        aMap[i,1] = (masterNodeNum-1)*6 + activeDof[i]
    end

    #determine global active DOFs associated with slave node
    aMap2 =zeros(length(slaveActiveDof))
    for i=1:length(slaveActiveDof)
        aMap2[i] = (slaveNodeNum-1)*6 + slaveActiveDof[i]
    end

    if (!isempty(aMap2)) #create overall map of active DOFs associated with this joint
        aDOF = [aMap;aMap2]
    else
        aDOF = aMap
    end
    return Tda,dDOF,aDOF
end

"""
Internal, determines the local master DOF associated with a local slave DOF.
"""
function determineActiveDofsFromSlaveNode(slaveDof,numDofPerNode)
    # Get size
    count = 0
    for i=1:numDofPerNode #loop over number of DOF per node
        if !in(i,slaveDof) #if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node
            #         slaveNodeActiveDof(count) = i;
            count = count + 1;

        end
    end

    if count>0
        slaveNodeActiveDof = zeros(count)
        count = 0
        for i=1:numDofPerNode #loop over number of DOF per node
            if !in(i,slaveDof) #if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node
                count = count + 1
                slaveNodeActiveDof[count] = i
            end
        end
    else
        slaveNodeActiveDof = []
    end

    return slaveNodeActiveDof
end

"""
Internal, gets the total number of DOFs in the model, active number of DOFs in the model, and a list of slave DOFs that will be eliminated by joint constraints.
"""
function extractdaInfo(joint,numNodes,numDofPerNode)

    adNumDof = numNodes*numDofPerNode; #total number of DOFs (active and dependent)

    numJoints=size(joint)[1] #get number of joints

    #Get Count
    dependentCount = 0
    count = 1
    for i=1:numJoints #loop over number of joints
        if (joint[i,4]==0 || joint[i,4]==5) #for a "fixed/weld" joint type or rigid bar constraint
            for j=1:6
                count = count + 1;
            end

        end

        if (joint[i,4]==1) #for a "pinned" joint type
            for j=1:3
                count = count + 1;
            end
        end

        if (joint[i,4]==2) #for a single axis hinge joint along a local "2" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

        if (joint[i,4]==3) #for a single axis hinge =joint along a local "1" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

        if (joint[i,4]==4) #for a single axis hinge joint along a local "3" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

    end

    slaveDof = zeros(count)
    count = 1
    for i=1:numJoints #loop over number of joints
        if (joint[i,4]==0 || joint[i,4]==5) #for a "fixed/weld" joint type or rigid bar constraint
            con = [1 2 3 4 5 6] #all DOFs of a slave node are constrained
            dependentCount = dependentCount + 6 #increment number of dependent DOFs
            for j=1:6
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end

        end

        if (joint[i,4]==1) #for a "pinned" joint type
            con = [1 2 3] #only translational (first 3) DOFs of a slave node are  constrained
            dependentCount = dependentCount + 3 #increment number of dependent DOFs
            for j=1:3
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==2) #for a single axis hinge joint along a local "2" axis of a joint
            if ((abs(abs(joint[i,7])-90))<1.0e-3 || (abs(abs(joint[i,7])-270))<1.0e-3)
                con=[1 2 3 5 6]
            else
                con=[1 2 3 4 6] #all but 5th DOF of a  slave node are constrained
            end
            dependentCount = dependentCount + 5 #increment number of dependent DOFs
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==3) #for a single axis hinge =joint along a local "1" axis of a joint
            if ((abs(abs(joint[i,8])-90))<1.0e-3 || (abs(abs(joint[i,8])-270))<1.0e-3)
                con = [1 2 3 4 5]
            elseif ((abs(abs(joint[i,7])-90))<1.0e-3 || (abs(abs(joint[i,7])-270))<1.0e-3)
                con = [1 2 3 4 6]
            else
                con = [1 2 3 5 6] #all but the 4th DOF of a slave node are constrained
            end
            dependentCount = dependentCount + 5 #increment number of dependent DOFs
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==4) #for a single axis hinge joint along a local "3" axis of a joint
            if ((abs(abs(joint[i,8])-90))<1.0e-3 || (abs(abs(joint[i,8])-270))<1.0e-3)
                con = [1 2 3 5 6]
                if ((abs(abs(joint(i,7))-90))<1.0e-3 || (abs(abs(joint(i,7))-270))<1.0e-3)
                    con = [1 2 3 4 6]
                end
            else
                con = [1 2 3 4 5]
            end
            dependentCount = dependentCount + 5 #all but the 6th DOF of a slave node are constrained
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

    end

    aNumDof = adNumDof - dependentCount #calculate number of active DOFs in the model

    return adNumDof,aNumDof,slaveDof

end

"""

   createJointTransform(joint,numNodes,numDofPerNode)

Internal, calculates the JointTransform of a structural system.

#Input
* `joint`:         object containing joint data
* `numModes`:      number of nodes in mesh
* `numDofPerNode`: number of degrees of freedom per node

#Output
* `jointTransform`: joint transformation matrix
* `reducedDOF`:     map of original DOF numbering to reduced DOF numbering
"""
function createJointTransform(joint,numNodes,numDofPerNode)

    numJoints=size(joint)[1]  #get number of joints in model

    #extract number of active DOFs, number of dependent DOFs, slave DOF numbers
    adNumDof,aNumDof,slaveDof = extractdaInfo(joint,numNodes,numDofPerNode)

    #initialize joint transformation matrix
    jointTransform = zeros(adNumDof,aNumDof)

    #form reduced DOF vector which maps original DOF numbering to reduced DOF
    #numbering

    #Get Count
    count = 1
    for i=1:numNodes*numDofPerNode #loop over total number of DOFs in model
        if !in(i,slaveDof)
            #         reducedDOF(count) = i #if DOF is NOT a slave DOF include it in reducedDOF
            count = count + 1
        end
    end

    reducedDOF = zeros(Int,count-1)
    count = 1
    for i=1:numNodes*numDofPerNode #loop over total number of DOFs in model
        if !in(i,slaveDof)
            reducedDOF[count] = i #if DOF is NOT a slave DOF include it in reducedDOF
            count = count + 1
        end
    end

    #create identity portion of transformation matrix (This is done by Craig,
    #but here the original DOF ordering is retained
    for i=1:aNumDof #loop over number of active DOFs
        jointTransform[reducedDOF[i],i] = 1.0 #mapping of active DOFs in full DOF list to reduced DOF list
    end

    #impose Tda portion of identity matrix and map to appropriate locations

    for i=1:numJoints # loop of number of joints in the model
        jointType = joint[i,4] #get joint type
        slaveNodeNum = joint[i,3] #get slave node number associated with joint
        masterNodeNum = joint[i,2] #get master node number associated with joint
        psi = joint[i,7] #get psi orientation angle associated with joint
        theta = joint[i,8] #get theta orientation angle associated with joint

        #Tda is a local transform between dependent and active DOFs for nodes
        #associated with a particular joint, dDOF is a listing of dependent
        #global DOFs associated with this joint, aDOF is a listing of
        #active global DOFs associated with this joint.
        Tda,dDOF,aDOF = createTda(jointType,slaveNodeNum,masterNodeNum,psi,theta,joint[i,:])

        for m=1:length(aDOF) #loop over global active DOFs associated with joint
            for k = 1:length(dDOF) #loop over global dependent DOFs associated with joint
                entry=findall(x->x==aDOF[m],reducedDOF)[1]  #determine reduced DOF associated with active DOF from original DOF listing
                jointTransform[dDOF[k],entry] = Tda[k,m]  #map local joint transformation matrix (Tda) to entries in global transformation matrix (jointTransform)
            end
        end
    end

    return jointTransform, reducedDOF
end

"""
Internal, searches over all DOFs in a structural model and determines and returns "dofVector" containing only unconstrained DOFs
"""
function calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained)

    #loop over all DOFs in the model checking if constrained by BC or not
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                #             dofVector(index) = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    dofVector = zeros(Int,index)
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                dofVector[index] = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    return dofVector
end

"""
Internal, creates a map of unconstrained DOFs between a full listing and reduced listing (after constraints have been applied)
"""
function constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,numpBC,pBC,isConstrained)

    bcdoflist=zeros(Int, numpBC)

    #create a listing of constrained DOFs from boundary condition file
    for i=1:numpBC
        bcnodenum = pBC[i,1]
        bcdofnum = pBC[i,2]
        bcdoflist[i] = (bcnodenum-1)*numDofPerNode + bcdofnum
    end

    dofList = calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained) #calculate a reduced (unconstrained) DOF vector

    redVectorMap = zeros(numReducedDof)

    for i=1:numReducedDof

        if (i in bcdoflist)              #creates a map of unconstrained reduced DOFs
            redVectorMap[i] = -1.0
        else
            index = findall(x->x==i,dofList)[1]
            redVectorMap[i] = index
        end

    end
    return redVectorMap
end


"""

   calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)

Internal, creates a boundary condition map between full and reduced dof listing as a result of constraints.

#Input
* `numpBC`            number of boundary conditions
* `pBC`               array of boundary  condition data
* `numDofPerNode`     number of degrees of freedom per node
* `reducedDofList`    array of reduced DOF numbering

#Output
* `elStorage`         map for boundary conditions between full and reduced dof list
"""
function calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)

    constrainedDof = zeros(numpBC)
    for i=1:numpBC
        constrainedDof[i] = (pBC[i,1]-1)*numDofPerNode + pBC[i,2]  #creates an array of constrained DOFs
    end
    constrainedDof = sort(constrainedDof)

    reducedDOFCount = length(reducedDofList)

    bcMap = zeros(reducedDOFCount)
    index = 1
    for i=1:reducedDOFCount
        if reducedDofList[i] in constrainedDof  #searches reduced DOF for constrained DOFs
            bcMap[i] = -1
        else
            bcMap[i] = index
            index = index + 1
        end
    end

    return bcMap

end

"""

   calculateShapeFunctions(elementOrder,xi,x)

This function calculates the Lagrange shape function, shape
function derivative, and Jacobian to map between the local element
domain and physical length of the element. The shape function
derivative is defined with respect to the physical length domain. The
shape functions may be linear or quadratic in order.

#Input
* `elementOrder` order of element: 1 linear, 2 quadratic
* `xi`           guass point values to evaluate shape functions at
* `x`            nodal coordinates in physical length domain

#Output
* `N`            shape function value at specified gauss points
* `p_N_x`        shape function derivative w.r.t physical length domain at specified gauss points
* `Jac`          Jacobian for mat between local element domain and physical length domain.
"""
function calculateShapeFunctions(elementOrder,xi,x)

    # N shape function
    # p_N_xi partial derivative of shape function w.r.t. xi

    #Linear interpolation functions
    N = zeros(elementOrder+1)
    p_N_xi = zeros(elementOrder+1)
    if elementOrder == 1
        N[1] = 0.5*(1.0 - xi)
        N[2] = 0.5*(1.0 + xi)

        p_N_xi[1] = -0.5
        p_N_xi[2] = 0.5
    end

    #Quadratic interpolation functions
    if elementOrder == 2
        N[1] = 0.5*(xi-1.0)*xi
        N[2] = 1.0-xi^2
        N[3] = 0.5*(xi+1.0)*xi

        p_N_xi[1] = xi - 0.5
        p_N_xi[2] = -2.0*xi
        p_N_xi[3] = xi + 0.5
    end

    numNodesPerEl = length(N)
    Jac=0.0
    for i=1:numNodesPerEl
        Jac = Jac + p_N_xi[i]*x[i]
    end

    p_N_x = zeros(numNodesPerEl)
    for i=1:numNodesPerEl
        p_N_x[i] = p_N_xi[i]/Jac
    end
    return N,p_N_x,Jac
end

"""
Internal, linear interpolation
"""
function interpolateVal(valNode,N)
    valGP = 0.0
    for i=1:length(N)
        valGP = valGP + N[i]*valNode[i]
    end
    return valGP
end

"""
Internal, forms total force vector and transform to desired DOF mapping
"""
function mapVector(Ftemp)

    a=length(Ftemp)
    Fel=zeros(a)

    # #declare map
    map = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]

    for i=1:a
        I=map[i]
        Fel[I] = Ftemp[i]
    end
    return Fel
end

"""
Internal, general routine to calculate an element matrix
"""
function calculateElement1!(EA,integrationFactor,N1,N2,K)
    len1 = length(N1)
    len2 = length(N2)
    for i=1:len1
        for j=1:len2
            K[i,j] = K[i,j] + EA*N1[i]*N2[j]*integrationFactor
        end
    end
end

"""
Internal, general routine to calculate an element vector
"""
function calculateVec1!(f,integrationFactor,N,F)
    len=length(N)
    for i=1:len
        F[i] = F[i] + f*N[i]*integrationFactor
    end
end

"""

    getGP(numGP)

Internal, defines gauss point coordinates in a local element frame
and the associated weights for Gaussian quadrature numerical
integration.

#Input
* `numGP`:  number of quad points used for integration

#Output
* `xi`:     list of quad point coordinates in local element frame
* `weight`: associated weights for quad point coordinate
"""
function getGP(numGP)

    #define Gauss integration points
    xi = zeros(numGP)
    weight = zeros(numGP)
    if (numGP == 1)
        xi[1] = 0
        weight[1] = 2.0
    elseif (numGP == 2)
        xi[1] = -sqrt(1/3)
        xi[2] = sqrt(1/3)
        weight[1] = 1.0
        weight[2] = 1.0
    elseif (numGP == 3)
        xi[1] = -sqrt(3.0/5.0)
        xi[2] = 0.0
        xi[3] = sqrt(3.0/5.0)
        weight[1] = 5.0/9.0
        weight[2] = 8.0/9.0
        weight[3] = 5.0/9.0
    elseif (numGP == 4)
        xi[1] = sqrt((3.0-2*sqrt(6.0/5.0))/7.0)
        xi[2] = -sqrt((3.0-2*sqrt(6.0/5.0))/7.0)
        xi[3] = sqrt((3.0+2*sqrt(6.0/5.0))/7.0)
        xi[4] = -sqrt((3.0+2*sqrt(6.0/5.0))/7.0)

        weight[1] = (18+sqrt(30))/36.0
        weight[2] = (18+sqrt(30))/36.0
        weight[3] = (18-sqrt(30))/36.0
        weight[4] = (18-sqrt(30))/36.0
    end
    return xi, weight
end

"""
Internal, calculates element mass properties.
"""
function calculateElementMass(rhoA,rhoIyy,rhoIzz,rhoIyz,rhoJ,ycm,zcm,x,y,z,integrationFactor,M,Itens,xm)
    M += rhoA*integrationFactor
    y += ycm
    z += zcm

    #Total MOI = parallel axis theorem + local MOI
    Itens .+= rhoA.*integrationFactor.*[(y^2+z^2)  -x*y  -x*z
    -x*y  (x^2+z^2) -y*z
    -x*z -y*z (x^2+y^2)] + integrationFactor.*[rhoJ 0 0
    0 rhoIyy rhoIyz
    0 rhoIyz rhoIzz]

    xm[1] += x*rhoA*integrationFactor
    xm[2] += y*rhoA*integrationFactor
    xm[3] += z*rhoA*integrationFactor
    return M,Itens,xm
end

"""

    ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,nodalStiffnessTerms,nodalLoads)

Compiles concentrated mass, stiffness, and load
associated with a node from both ndl and joint files. The mod*
variables are passed back with these terms removed to prevent
duplicate application of shared nodal terms between elements

#Input
* `conn`                connectivity list for element
* `joint`               joint array for nodal terms
* `nodalMassTerms`      listing of concentrated nodal mass terms
* `nodalStiffnessTerms` listing of concentrated nodal stiffness terms
* `nodalLoads`          listing of concentrated nodal loads terms


#Output
* `mass`                array of concentrated mass associated with element
* `stiff`               array of concentrated stiffness associated with element
* `load`                array of concentrated loads associated with element
* `modJoint`            modified joint object removing nodal terms that have/will be applied to the element calculations
* `modNodalMassTerms`   modified nodal mass object removing nodal terms that have/will be applied to the element calculations
* `modalStiffnessTerms` modified nodal stiffness object removing nodal terms that have/will be applied to the element calculations
* `modNodalLoads`       modified nodal loads object removing nodal terms that have/will be applied to the element calculations

"""
function ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,nodalStiffnessTerms,nodalLoads)

    node1 = conn[1] #define node #1 and node #2
    node2 = conn[2]

    mass1=0  #initialize concentrated mass amd moi for nodes
    mass2=0
    moix1=0
    moiy1=0
    moiz1=0
    moix2=0
    moiy2=0
    moiz2=0

    stiff1x=0 #initialize concentrated stifness for nodes
    stiff2x=0
    stiff1y=0
    stiff2y=0
    stiff1z=0
    stiff2z=0
    stiff1mx=0
    stiff2mx=0
    stiff1my=0
    stiff2my=0
    stiff1mz=0
    stiff2mz=0

    f1x = 0   #initialize concentrated loads/moments
    f2x = 0
    f1y = 0
    f2y = 0
    f1z = 0
    f2z = 0
    m1x =0
    m2x =0
    m1y =0
    m2y =0
    m1z =0
    m2z =0

    modJoint = joint                         #create copies of joint, and nodal mass, stiffness, loads arrays
    modNodalMassTerms = deepcopy(nodalMassTerms)
    modNodalStiffnessTerms = deepcopy(nodalStiffnessTerms)
    modNodalLoads = deepcopy(nodalLoads)

    numJoints,_=size(joint)    #get number of joints in model

    if numJoints > 0
        node1flag=joint[:,2].==node1  #see if nodes are associated with a joint constraint as a master node
        node2flag=joint[:,2].==node2
    else
        node1flag = false
        node2flag = false
        mass1 = 0.0
        mass2 = 0.0
    end

    for i=1:numJoints           #if nodes are associated with joint constraint, use (if any) mass and stiffness specification from the joint file
        if node1flag[i]==1
            mass1 = mass1+joint[i,5]
            #             stiff1x = stiff1x + joint[i,6]
            #             stiff1y = stiff1y + joint[i,6]
            #             stiff1z = stiff1z + joint[i,6]
            #             stiff1mx = stiff1mx + joint[i,6]
            #             stiff1my = stiff1my + joint[i,6]
            #             stiff1mz = stiff1mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
        if node2flag[i]==1
            mass2 = mass2+joint[i,5]
            #             stiff2x = stiff2x + joint[i,6]
            #             stiff2y = stiff2y + joint[i,6]
            #             stiff2z = stiff2z + joint[i,6]
            #             stiff2mx = stiff2mx + joint[i,6]
            #             stiff2my = stiff2my + joint[i,6]
            #             stiff2mz = stiff2mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
    end

    #apply concentrated mass/stiffness from NDL file

    for i=1:length(nodalMassTerms)   #if node is specified in nodal mass terms file add to mass properties for this node
        node1flagM=nodalMassTerms[i].nodeNum.==node1
        node2flagM=nodalMassTerms[i].nodeNum.==node2
        if node1flagM==1
            if nodalMassTerms[i].dof == 1
                mass1 = mass1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 4
                moix1 = moix1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy1 = moiy1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz1 = moiz1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
        if node2flagM==1
            if nodalMassTerms[i].dof == 1
                mass2 = mass2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 4
                moix2 = moix2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy2 = moiy2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz2 = moiz2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
    end



    for i=1:length(nodalStiffnessTerms)     #if node is specified in nodal stiffness terms file add to stiffness properties for this node
        node1flagK=nodalStiffnessTerms[i].nodeNum.==node1
        node2flagK=nodalStiffnessTerms[i].nodeNum.==node2
        if node1flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff1x = stiff1x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff1y = stiff1y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff1z = stiff1z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff1mx = stiff1mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff1my = stiff1my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff1mz = stiff1mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
        if node2flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff2x = stiff2x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff2y = stiff2y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff2z = stiff2z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff2mx = stiff2mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff2my = stiff2my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff2mz = stiff2mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
    end

    for i=1:length(nodalLoads)  #if node is specified in nodal forces terms file add to concentrated force for this node
        node1flagF=nodalLoads[i].nodeNum.==node1
        node2flagF=nodalLoads[i].nodeNum.==node2
        if node1flagF==1
            if nodalLoads[i].dof==1
                f1x = f1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f1y = f1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f1z = f1z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m1x = m1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m1y = m1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m1z = m1z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
        if node2flagF==1
            if nodalLoads[i].dof==1
                f2x = f2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f2y = f2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f2z = f2z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m2x = m2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m2y = m2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m2z = m2z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
    end


    #compile nodal concentrated terms into mass, stiffness, and load arrays
    mass = [mass1 mass2
    moix1 moix2
    moiy1 moiy2
    moiz1 moiz2]

    stiff = [stiff1x stiff2x
    stiff1y stiff2y
    stiff1z stiff2z
    stiff1mx stiff2mx
    stiff1my stiff2my
    stiff1mz stiff2mz]

    load = [f1x f2x
    f1y f2y
    f1z f2z
    m1x m2x
    m1y m2y
    m1z m2z]

    return mass,stiff,load,modJoint,modNodalMassTerms,modNodalStiffnessTerms,modNodalLoads

end

"""

    assembly(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg)

Internal, assembles the element matrix and load vector into the
global system of equations

#Input
* `Ke`:            element matrix
* `Fe`:            element vector
* `conn`:          element connectivity
* `numNodesPerEl`: number of nodes per element
* `numDofPerNode`: number of degrees of freedom per node
* `Kg`:            global system matrix
* `Fg`:            global load vector

#Output
* `Kg`:            global system matrix with assembled element
* `Fg`:            global load vector with assembled element

"""
function assembly!(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg)

    count = 1
    dofList = zeros(Int,numNodesPerEl*numDOFPerNode)
    for i=1:numNodesPerEl
        for j=1:numDOFPerNode
            dofList[count] = (conn[i]-1)*numDOFPerNode + j
            count = count + 1
        end
    end

    numDOFPerEl = length(dofList)
    #Assemble element i into global system
    for j=1:numDOFPerEl
        J = dofList[j]
        Fg[J] = Fg[J] + Fe[j]
        for m=1:numDOFPerEl
            M = dofList[m]
            Kg[J,M] = Kg[J,M] + Ke[j,m]
        end
    end
end

function assembly(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg)

    count = 1
    dofList = zeros(Int,numNodesPerEl*numDOFPerNode)
    for i=1:numNodesPerEl
        for j=1:numDOFPerNode
            dofList[count] = (conn[i]-1)*numDOFPerNode + j
            count = count + 1
        end
    end

    numDOFPerEl = length(dofList)
    #Assemble element i into global system
    for j=1:numDOFPerEl
        J = dofList[j]
        Fg[J] = Fg[J] + Fe[j]
        for m=1:numDOFPerEl
            M = dofList[m]
            Kg[J,M] = Kg[J,M] + Ke[j,m]
        end
    end
    return Kg,Fg
end

"""

    applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode)

Internal, applies boundary conditions to the stiffness matrix and
load vector for a static analysis.

#Input
* `Kg`            assembled global stiffness matrix
* `Fg`            assembled global load vector
* `BC`            struct of boundary condition information
* `u`             global displacement vector
* `iterationType` for nonlinear analysis, not used in BLAST
* `numDofPerNode` number of degrees of freedom per node

#Output
* `Kg`            global stiffness matrix with boundary conditions
* `Fg`            global load vector with boundary condition

"""
function applyBC(Kg,Fg,BC,numDofPerNode)

    numEq=size(Kg)[1]

    #APPLY BCs FOR PRIMARY VARIABLE

    if (BC.numpBC > 0)
        pBC = BC.pBC
        numpBC = size(pBC)[1]
        eqidx = findall(x->x==-1,BC.map)
        for i=1:numpBC
            nodeNumber = Int(pBC[i,1])
            dofNumber = Int(pBC[i,2])
            specVal = pBC[i,3]

            eqNumber = eqidx[i]#(nodeNumber-1)*numDofPerNode + dofNumber

            for j=1:numEq
                Kg[eqNumber,j] = 0.0
                Fg[j] = Fg[j] - Kg[j,eqNumber]*specVal
                Kg[j,eqNumber] = 0.0
            end
            Fg[eqNumber] = specVal
            Kg[eqNumber,eqNumber] = 1.0
        end
    end

    #APPLY BCs FOR SECONDARY VARIABLE

    # if (BC.numsBC > 0) # This does not appear to be used
    #     sBC = BC.sBC
    #     [numsBC,~] = size(sBC)
    #
    #     for i=1:numsBC
    #         nodeNumber = sBC[i,1]
    #         dofNumber = sBC[i,2]
    #         specVal =  sBC[i,3]
    #
    #         eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber
    #
    #         Fg(eqNumber) = Fg(eqNumber) + specVal
    #
    #     end
    # end
    return Kg,Fg
end

"""

    calculateReactionForceAtNode(nodeNum,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H,countedNodes)

Internal, calculates the reaction force at a node by post
processing all element associated with a node through connectivity or
joint constraints.

#Input
* `nodeNum`:    node number joint constraints are desired at
* `model`:      object containing model data
* `mesh`:       object containing mesh data
* `elStorage`:  object containing stored element data
* `el`:         object containing element data
* `timeInt`:    object containing time integration parameters
* `dispData`:   object containing displacement data
* `displ_iter`: converged displacement solution
* `rbData`:     vector containing rigid body displacement, velocity, and acceleration
* `Omega`:      rotor speed (Hz)
* `OmegaDot`:   rotor acceleratin (Hz)
* `CN2H`:       transformation matrix from inertial frame to hub frame
* `countedNodes`:  prevents nodal terms from being double counted


#Output
* `cummulativeForce`:  vector containing reaction force at nodeNum
"""
function calculateReactionForceAtNode(nodeNum,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H,countedNodes)

    conn = mesh.conn  #get connectivity list
    numDofPerNode = 6

    cummulativeForce = zeros(numDofPerNode) #initialize force at node

    #find elements associated with nodeNum due to mesh connectivity or
    #constraints
    elList,elListLocalNodeNumbers = findElementsAssociatedWithNodeNumber(nodeNum,conn,model.joint)

    #process elements for nodal reaction forces and compile to find total
    #reaction force at specified node
    for i=1:length(elList)
        Fpp = elementPostProcess(elList[i],model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H,countedNodes)
        localNode = elListLocalNodeNumbers[i]
        cummulativeForce = cummulativeForce + Fpp[(localNode-1)*numDofPerNode+1:(localNode-1)*numDofPerNode+6]
    end
    return cummulativeForce
end

"""
Internal calculates element strains
"""
function calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ,nlflag)

    elStrain = Array{ElStrain, 1}(undef, numEl)

    for i=1:numEl
        #Calculate Ke and Fe for element i
        index = 1
        # elementOrder = elementOrder
        # nlOn = nlflag
        xloc = [0.0 el.elLen[i]]
        sectionProps = el.props[i]
        sweepAngle = el.psi[i]
        coneAngle = el.theta[i]
        rollAngle = el.roll[i]
        aeroSweepAngle = 0.0

        eldisp = zeros(1,numNodesPerEl*numDOFPerNode)
        for j=1:numNodesPerEl       #define element coordinates and displacements associated with element
            for k=1:numDOFPerNode
                eldisp[index] = displ[(conn[i,j]-1)*numDOFPerNode + k] #TODO: add in logic that if it is on a constrained/joint node that it uses the displacement from the prior?
                index = index + 1
            end
        end

        # disp = eldisp
        elStrain[i] = calculateTimoshenkoElementStrain(elementOrder,nlflag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,eldisp)

    end
    return elStrain
end

"""

    findElementsAssociatedWithNodeNumber(nodeNum,conn,jointData)

Internal, finds elements associated with a node number through
mesh connectivity or joint constraints

#Input
* `nodeNum`    node number joint constraints are desired at
* `conn`       object containing mesh connectivity
* `jointData`  object containing joint information


#Output
* `elList`     array containing a list of element numbers associated with nodeNum
* `localNode`  array containing the local node number that correspond to nodeNum in the list of associated elements
"""
function findElementsAssociatedWithNodeNumber(nodeNum,conn,jointData)

    #search joint constraints
    index = 1
    numEl = size(conn)[1] #get number of elements in mesh
    elList = zeros(Int,0)
    localNode = zeros(Int,0)
    if !isempty(jointData)
        #first see if specified node is a slave node in a joint constraint
        # keep this here for future translation from matlab: res2 = find(ismember(jointData(:,3),nodeNum)) #search joint data slave nodes for node number
        #if it is, change it to the corresponding master node
        if !isempty(findall(x->x==nodeNum,jointData[:,3]))
            nodeNum = jointData[end,2]
            if length(jointData)>1
                error("Incorrect Joint Data and nodeNum, too many joints")
            end
        end

        res1 = findall(x->x==nodeNum,jointData[:,2]) #search joint data master nodes for node number
        if !isempty(res1)
            jointNodeNumbers = jointData[res1,3]

            for j=1:length(jointNodeNumbers) #loop over joints
                for i=1:numEl
                    localNodeNumber = findall(x->x==jointNodeNumbers[j],conn[i,:]) #finds indices of nodeNum in connectivity of element i #finds the local node number of element i that corresponds to nodeNum
                    if !isempty(localNodeNumber) #assigns to an elementList and localNode list
                        elList = vcat(elList,i)
                        localNode = vcat(localNode,localNodeNumber[1])
                        index = index + 1
                    end
                end
            end
        end
    # else
    #     error("empty jointData")
    end


    for i=1:numEl #loop over elements
        localNodeNumber = findall(x->x==nodeNum,conn[i,:]) #finds indices of nodeNum in connectivity of element i #finds the local node number of element i that corresponds to nodeNum
        if !isempty(localNodeNumber) #assigns to an elementList and localNode list
            elList = vcat(elList,i)
            localNode = vcat(localNode,localNodeNumber[1])
            index = index + 1
        end
    end

    return elList,localNode

end

"""
Internal, gets the concentrated terms without double counting
"""
function getElementConcTerms!(Kconc, Mconc, Cconc, Fconc, elNodes, numDOFPerNode, appliedNodes)

    if elNodes[1] in appliedNodes
        elKconc1 = zeros(numDOFPerNode, numDOFPerNode)
        elMconc1 = zeros(numDOFPerNode, numDOFPerNode)
        elCconc1 = zeros(numDOFPerNode, numDOFPerNode)
        elFconc1 = zeros(numDOFPerNode)
    else
        elKconc1 = Kconc[(elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode, (elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode]
        elMconc1 = Mconc[(elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode, (elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode]
        elCconc1 = Cconc[(elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode, (elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode]
        elFconc1 = Fconc[(elNodes[1]-1)*numDOFPerNode+1:elNodes[1]*numDOFPerNode]
        append!(appliedNodes, elNodes[1])
    end

    if elNodes[2] in appliedNodes
        elKconc2 = zeros(numDOFPerNode, numDOFPerNode)
        elMconc2 = zeros(numDOFPerNode, numDOFPerNode)
        elCconc2 = zeros(numDOFPerNode, numDOFPerNode)
        elFconc2 = zeros(numDOFPerNode)
    else
        elKconc2 = Kconc[(elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode, (elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode]
        elMconc2 = Mconc[(elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode, (elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode]
        elCconc2 = Cconc[(elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode, (elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode]
        elFconc2 = Fconc[(elNodes[2]-1)*numDOFPerNode+1:elNodes[2]*numDOFPerNode]
        append!(appliedNodes, elNodes[2])
    end

    elKconc = hcat(elKconc1, elKconc2)
    elMconc = hcat(elMconc1, elMconc2)
    elCconc = hcat(elCconc1, elCconc2)
    elFconc = hcat(elFconc1, elFconc2)

    return elKconc, elMconc, elCconc, elFconc, appliedNodes
end

"""

    writeOutput(freq,damp,phase1,phase2,imagComponentSign,fid)

Internal, writes an output file and or formats an output for modal analysis.

#Input
* `freq`:               array of modal frequencies
* `damp`:               array of modal damping ratios
* `phase1`:             array of in phase mode shapes
* `phase2`:             array of out of phase mode shapes
* `imagComponentSign`:  array of sign of imaginary components
* `fid`:                file identifier for output

#Output
* `freqSorted`:         array of sorted(by frequency) modal frequencies
* `dampSorted`:         array of sorted(by frequency) modal damping ratios
* `imagCompSignSorted`: array of sorted(by frequency) of imaginarycomponentSign array
* `U_x_0`: see ?Modal outputs
* `U_y_0`: see ?Modal outputs
* `U_z_0`: see ?Modal outputs
* `theta_x_0`: see ?Modal outputs
* `theta_y_0`: see ?Modal outputs
* `theta_z_0`: see ?Modal outputs
* `U_x_90`: see ?Modal outputs
* `U_y_90`: see ?Modal outputs
* `U_z_90`: see ?Modal outputs
* `theta_x_90`: see ?Modal outputs
* `theta_y_90`: see ?Modal outputs
* `theta_z_90`: see ?Modal outputs
"""
function ModalOutput(freq,damp,phase1,phase2,imagComponentSign,filename)

    if filename!="none"
        println("Saving Modal Output at: $filename")
        fid=open(filename,"w")
    end
    Nnodes = size(phase1[:,:,1])[1] #gets number of nodes for mode shape printing

    posIndex = 1
    Nmodes = length(freq)

    dampSorted = zeros(Nmodes)
    freqSorted = zeros(Nmodes)
    imagCompSignSorted = zeros(Nmodes)
    U_x_0 = zeros(Nnodes,Nmodes)
    U_y_0 = zeros(Nnodes,Nmodes)
    U_z_0 = zeros(Nnodes,Nmodes)
    theta_x_0 = zeros(Nnodes,Nmodes)
    theta_y_0 = zeros(Nnodes,Nmodes)
    theta_z_0 = zeros(Nnodes,Nmodes)
    U_x_90 = zeros(Nnodes,Nmodes)
    U_y_90 = zeros(Nnodes,Nmodes)
    U_z_90 = zeros(Nnodes,Nmodes)
    theta_x_90 = zeros(Nnodes,Nmodes)
    theta_y_90 = zeros(Nnodes,Nmodes)
    theta_z_90 = zeros(Nnodes,Nmodes)

    index = 1
    for i=posIndex:1:posIndex+(Nmodes-1) #prints mode frequency, damping and in/out of phase mode shapes
        if filename!="none"
            Printf.@printf(fid,"MODE # %0.0f \n\n",index)
            Printf.@printf(fid,"Frequency: %e: \n",freq[i])
            Printf.@printf(fid,"Damping %e: \n",damp[i])
            Printf.@printf(fid,"0 deg Mode Shape:\n")
            Printf.@printf(fid,"U_x          U_y          U_z          theta_x     theta_y     theta_z \n")
        end
        for j=1:Nnodes
            U_x_0[j,index] = phase1[j,1,i]
            U_y_0[j,index] = phase1[j,2,i]
            U_z_0[j,index] = phase1[j,3,i]
            theta_x_0[j,index] = phase1[j,4,i]
            theta_y_0[j,index] = phase1[j,5,i]
            theta_z_0[j,index] = phase1[j,6,i]
            if filename!="none"
                Printf.@printf(fid,"%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",U_x_0[j,index],U_y_0[j,index],U_z_0[j,index],theta_x_0[j,index],theta_y_0[j,index],theta_z_0[j,index])
                Printf.@printf(fid,"\n")
            end
        end

        if filename!="none"
            Printf.@printf(fid,"\n")

            Printf.@printf(fid,"90 deg Mode Shape:\n")
            Printf.@printf(fid,"U_x          U_y          U_z          theta_x     theta_y     theta_z \n")
        end

        for j=1:Nnodes
            U_x_90[j,index] = phase2[j,1,i]
            U_y_90[j,index] = phase2[j,2,i]
            U_z_90[j,index] = phase2[j,3,i]
            theta_x_90[j,index] = phase2[j,4,i]
            theta_y_90[j,index] = phase2[j,5,i]
            theta_z_90[j,index] = phase2[j,6,i]
            if filename!="none"
                Printf.@printf(fid,"%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",U_x_90[j,index],U_y_90[j,index],U_z_90[j,index],theta_x_90[j,index],theta_y_90[j,index],theta_z_90[j,index])
                Printf.@printf(fid,"\n")
            end
        end

        if (i<posIndex+(Nmodes-1)) && filename!="none"
            Printf.@printf(fid,"\n\n")
        end



        dampSorted[i] = damp[i]
        freqSorted[i] = freq[i]
        imagCompSignSorted[i] = imagComponentSign[i]

        index = index + 1
    end

    if filename!="none"
        close(fid)
    end

    return freqSorted,dampSorted,imagCompSignSorted,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
end

"""

   calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)

Takes in a global 6dof distributed force at two nodal points and returns the 6dof force in the element FOR

#Input
* `elementOrder::`:
* `x::Array{<:float}`: mesh x-position
* `xloc::Array{<:float}`: local x-position [0 elLength]
* `twist::Array{<:float}`: element twist angle (rad)
* `sweepAngle::Array{<:float}`: element sweep angle (rad)
* `coneAngle::Array{<:float}`: element cone angle (rad)
* `rollAngle::Array{<:float}`: element roll angle (rad)
* `extDistF2Node::Array{<:float}`: turbine Tangential force
* `extDistF3Node::Array{<:float}`: turbine Normal force
* `extDistF4Node::Array{<:float}`: turbine M25 moment

#Output
* `Fe::Array{float}`: 6x1 Force on element in element FOR
"""
function calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)

    numGP = 4   #number of gauss points for full integration
    #calculate quad points
    xi,weight = getGP(numGP)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    F1 = zeros(numNodesPerEl,1)
    F3 = zero(F1)
    F2 = zero(F1)
    F4 = zero(F1)
    F5 = zero(F1)
    F6 = zero(F1)

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg = rollAngle + 0.5*(twist[1] + twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,_,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = N
        N2 = N
        N3 = N
        N4 = N
        N5 = N
        N6 = N
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        extDistF1 = 0
        extDistF2 = interpolateVal(extDistF2Node,N2)
        extDistF3 = interpolateVal(extDistF3Node,N3)
        extDistF4 = interpolateVal(extDistF4Node,N4)
        extDistF5 = 0
        extDistF6 = 0

        #distributed/body force load calculations
        calculateVec1!(extDistF1,integrationFactor,N1,F1)
        calculateVec1!(extDistF2,integrationFactor,N2,F2)
        calculateVec1!(extDistF3,integrationFactor,N3,F3)
        calculateVec1!(extDistF4,integrationFactor,N4,F4)
        calculateVec1!(extDistF5,integrationFactor,N5,F5)
        calculateVec1!(extDistF6,integrationFactor,N6,F6)


    end #END OF INTEGRATION LOOP

    #compile element force vector
    Fe = mapVector([F1;F2;F3;F4;F5;F6])

    # transform matrices for sweep
    # Note,a negative sweep angle, will sweep away from the direction of
    # positive rotation
    lambdaTran = lambda'
    # lambdaTran = sparse(lambdaTran)
    Fe = lambdaTran*Fe

    return Fe

end

"""

    setInitialConditions(initCond,u,numDOFPerNode)

 sets initial conditions

#Input
* `initCond`: array containing initial conditions
                 initCond(i,1) node number for init cond i
                 initCond(i,2) local DOF number for init cond i
                 initCond(i,3) value for init cond i
* `u`: displacement vector for each dof
* `numDOFPerNode`: number of degrees of freedom per node

#Output
* `u`: displacement vector modified for initial conditions

"""
function setInitialConditions(initCond,u,numDOFPerNode)

    len=size(initCond) #get number of specified initial conditions
    #unspecified initial conditions are assumed to
    #be zero

    for i=1:len[1] #loop over initial conditions
        if (initCond[i,2]>numDOFPerNode) #error check
            error("setInitalConditios:: DOF greater than numDOFPerNode")
        end
        index = Int((initCond[i,1]-1)*numDOFPerNode + initCond[i,2]) #calculate global DOF number for initial condition
        u[index] = initCond[i,3] #specify initial condition using global DOF number
    end

    return u

end

"""

    makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)

Internal, usese the pBC matrix and calculates/stores boundary condition data

#Input
* `pBC`      See ?FEAModel.pBC
* `numNodes`        number of nodes in structural model
* `numDofPerNode`   number of degrees of freedom per node
* `reducedDOFList`  joint transformation matrix from reduced to full DOF list
* `jointTransform`  listing of reduced DOFs

#Output
* `BC:BC_struct` see ?BC_struct

"""
function makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)

    totalNumDof = numNodes*numDofPerNode

    numsBC = 0
    nummBC = 0

    #create a vector denoting constrained DOFs in the model (0 unconstrained, 1
    #constrained)

    #calculate constrained dof vector
    isConstrained = zeros(totalNumDof)
    constDof = (pBC[:,1].-1)*numDofPerNode + pBC[:,2]
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if ((i-1)*numDofPerNode + j in constDof)
                isConstrained[index] = 1
            end
            index = index + 1
        end
    end
    numpBC = length(pBC[:,1])

    map = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDOFList)
    numReducedDof = length(jointTransform[1,:])
    redVectorMap = constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,numpBC,pBC,isConstrained) #create a map between reduced and full DOF lists

    BC = BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    map,
    redVectorMap)

    return BC

end

"""

    applyConcentratedTerms(numNodes, numDOFPerNode; filename="none",data=[1 "M6" 1 1 0.0], jointData=[])

Internal, applies 6x6 concentrated nodal terms from user input.

#Input
* `filename`:   string containing nodal terms filename
* `data`:       Nx5 or Nx4 array matching general [1 "M6" 1 1 0.0] or diagonal only [1 "M" 1 0.0] (node, Type, dof, value) where type is M,C,K, or F

#Output
* `nodalTerms::NodalTerms`: see ?NodalTerms object containing concentrated nodal data
"""
function applyConcentratedTerms(numNodes, numDOFPerNode; filename="none",data=[], jointData=[])

    if filename!="none" && !occursin("[",filename)
        data = DelimitedFiles.readdlm(filename,' ',skipstart = 0)
    end

    concStiff = zeros(numNodes*numDOFPerNode, numNodes*numDOFPerNode)
    concMass = zeros(numNodes*numDOFPerNode, numNodes*numDOFPerNode)
    concDamp = zeros(numNodes*numDOFPerNode, numNodes*numDOFPerNode)
    concLoad = zeros(numNodes*numDOFPerNode)
    if !isempty(data)
        for i_data = 1:length(data[:,1])
            if data[i_data,2][1] == 'M'
                nodeNum = data[i_data,1]
                dof1 = data[i_data,3]
                if length(data[1,:])==5 #If 6x6 general method
                    dof2 = data[i_data,4]
                elseif length(data[1,:])==4
                    @warn "Only one degree of freedom given for concentrated mass; applying at diagonal."
                    dof2 = data[i_data,3]
                else
                    error("Wrong number of inputs in the concentrated mass term.")
                end
                val = data[i_data,end]
                gdof1 = (nodeNum-1)*6 + dof1
                gdof2 = (nodeNum-1)*6 + dof2
                concMass[gdof1,gdof2] += val

            elseif data[i_data,2][1] == 'K'
                nodeNum = data[i_data,1]
                dof1 = data[i_data,3]
                if length(data[1,:])==5 #If 6x6 general method
                    dof2 = data[i_data,4]
                elseif length(data[1,:])==4
                    @warn "Only one degree of freedom given for concentrated stiffness; applying at diagonal."
                    dof2 = data[i_data,3]
                else
                    error("Wrong number of inputs in the concentrated stiffness term.")
                end
                val = data[i_data,end]
                gdof1 = (nodeNum-1)*6 + dof1
                gdof2 = (nodeNum-1)*6 + dof2
                concStiff[gdof1,gdof2] += val

            elseif data[i_data,2][1] == 'C'
                nodeNum = data[i_data,1]
                dof1 = data[i_data,3]
                if length(data[1,:])==5 #If 6x6 general method
                    dof2 = data[i_data,4]
                elseif length(data[1,:])==4
                    @warn "Only one degree of freedom given for concentrated damping; applying at diagonal."
                    dof2 = data[i_data,3]
                else
                    error("Wrong number of inputs in the concentrated damping term.")
                end
                val = data[i_data,end]
                gdof1 = (nodeNum-1)*6 + dof1
                gdof2 = (nodeNum-1)*6 + dof2
                concDamp[gdof1,gdof2] += val

            elseif data[i_data,2][1] == 'F'
                nodeNum = data[i_data,1]
                if length(data[1,:])==5
                    @warn "Concentrated loads are one-dimensional, second degree of freedom input is ignored."
                elseif length(data[1,:])!=4
                    error("Wrong number of inputs in the concentrated load term.")
                end
                dof = data[i_data,3]
                val = data[i_data,end]
                gdof = (nodeNum-1)*numDOFPerNode + dof
                concLoad[gdof] += val
            else
                error("Unknown Nodal Data Type")
            end
            if size(jointData)[1] > 0 # if there are joints in the model, apply joint mass to the concentrated mass for its master node
                for i=1:size(jointData)[1]
                    nodeNum = Int.(jointData[i,2])
                    massdofs = [1 1; 2 2; 3 3] # the first three diagonal terms of the mass matrix are the pure mass at the node
                    val = jointData[i,5]
                    gdofs = Int.((nodeNum-1)*6 .+ massdofs)

                    for j=1:size(gdofs)[1]
                        concMass[gdofs[j,:]] .+= val
                    end
                end
            end
        end
    end

    #store concentrated nodal term data in nodalTerms object
    return NodalTerms(concLoad,concStiff,concMass,concDamp)

end
