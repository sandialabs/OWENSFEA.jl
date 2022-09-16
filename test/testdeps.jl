import Statistics
# mesh for a single continuous beam
# Rotating beam case from Finite Element Solution of Nonlinear Intrinsic Equations for Curved Composite Beams by Hodges, Shang, and Cesnic
# ______
#       \
function mesh_beam(;L1 = 31.5, #first section of beam length
    L2 = 6.0, #second section of beam length
    Nelem1 = 13,
    Nelem2 = 3,
    angleD = 45.0, # angle of second section of beam relative to first (0 is straight)
    zeroOffset = 2.5,
    vertical=true)

    angle = angleD/360*2*pi

    # First Section
    mesh_x1 = collect(LinRange(zeroOffset,zeroOffset+L1,Nelem1+1))
    mesh_y1 = zero(mesh_x1)
    mesh_z1 = zero(mesh_x1)

    # intra-beam connectivity
    conn1 = zeros(length(mesh_z1)-1,2)
    conn1[:,1] = collect(1:length(mesh_z1)-1)
    conn1[:,2] = collect(2:length(mesh_z1))

    # Angled Section
    mesh_x2 = zeroOffset+L1.+cos(angle).*collect(LinRange(0.0,L2,Nelem2+1))
    mesh_y2 = -sin(angle).*collect(LinRange(0.0,L2,Nelem2+1))
    mesh_z2 = zero(mesh_x2)

    # intra-beam connectivity
    conn2 = zeros(length(mesh_z2)-1,2)
    conn2[:,1] = length(mesh_z1).+collect(1:length(mesh_z2)-1)
    conn2[:,2] = length(mesh_z1).+collect(2:length(mesh_z2))

    conn = [conn1;conn2]

    if vertical
        mesh_z = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_x = [mesh_z1;mesh_z2]
    else
        mesh_x = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_z = [mesh_z1;mesh_z2]
    end

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:end] .= 0 #Tower

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,2)

    meshSeg[1] = Nelem1
    meshSeg[2] = Nelem2

    # Not used for the beam case
    structuralSpanLocNorm = []
    structuralNodeNumbers = []
    structuralElNumbers = []
    # end

    mymesh = GyricFEA.Mesh(round.(Int,nodeNum),round.(Int,numEl),round.(Int,numNodes),mesh_x,mesh_y,mesh_z,round.(Int,elNum),round.(Int,conn),meshtype,meshSeg,structuralSpanLocNorm,round.(Int,structuralNodeNumbers),round.(Int,structuralElNumbers))

    ####################################
    ##----------Joint Matrix----------##
    ####################################

    #Connect L1 tip to L2 base
    jointconn = [length(mesh_x1) length(mesh_x1)+1]

    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh) #TODO: consider getting rid of ort struct for simplification since it isn't used hardly at all
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,2],ort.elNum) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2]-1,ort.elNum) #TODO: we get away with this since the elements are increasing and there are no two point objects (and this is only a problem for the top of the tower connecting to the blade tops)
        end
        if length(elnum_of_joint)==0
            elnum_of_joint = findall(x->x==jointconn[jnt,1]-1,ort.elNum)
        end
        Psi_d_joint[jnt] = ort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = ort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

    #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint) zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    return mymesh, ort, myjoint
end

function rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)
    #rigidBodyRotation rotates a vector through a rotation sequence
    #   [H1,H2,H3] = rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)
    #
    #   This function performs a coordinate transformation from a local
    #   body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
    #   rotation sequence
    #
    #      input:
    #      B1        = array containing body frame 1 coordinates of points to be
    #                  mapped to the hub frame
    #      B2        = array containing body frame 2 coordinates of points to be
    #                  mapped to the hub frame
    #      B3        = array containing body frame 3 coordinates of points to be
    #                  mapped to the hub frame
    #     AngleArray = Array of angles for Euler rotation sequence
    #     AxisArray  = Array of axes for Euler rotation sequence
    #
    #      output:
    #      H1        = array containg hub frame 1 coordinates of points mapped
    #                  to the hub frame from body frame
    #      H2        = array containg hub frame 2 coordinates of points mapped
    #                  to the hub frame from body frame
    #      H3        = array containg hub frame 3 coordinates of points mapped
    #                  to the hub frame from body frame

    #This function performs a coordinate transformation from a local
    #body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
    #rotation sequence

    #That is CHtoB = [M3(SweepAngle)][M2(Theta)][M3(Psi)];

    #calculate coordinate transformation matrix from element frame to
    #hub frame (CBtoH)
    dcm = createGeneralTransformationMatrix(AngleArray,AxisArray)
    C = dcm'

    #transform body coordinatized vector to be coordinatized in the hub
    #frame
    H1 = C[1,1].*B1 + C[1,2].* B2 + C[1,3].*B3
    H2 = C[2,1].*B1 + C[2,2].* B2 + C[2,3].*B3
    H3 = C[3,1].*B1 + C[3,2].* B2 + C[3,3].*B3

    return H1,H2,H3
end

function calculatePsiTheta(v)
    #calculatePsiTheta calculates the orienation of a single element
    #   [Psi,Theta] = calculatePsiTheta(v)
    #
    #   This function calculates the orientation of a single element. A local
    #   element frame is related to a hub frame through a transformation matrix
    #   CHtoE (transforming a vector from an element frame E to a global frame
    #   H) such that CHtoE = [M2(Theta)]*[M3(Psi)]. Here [M2( )] is a direction
    #   cosine matrix about a 2 axis and [M3( )] is a direction cosine matrix
    #   about a 3 axis.
    #
    #      input:
    #      v          = vector from node 1 to node 2 of an element
    #
    #      output:
    #      Psi        = "3" angle for element orientation (deg)
    #      Theta      = "2" angle for element orientation (deg)
    #                   see above for definition

    v = v./LinearAlgebra.norm(v) #normalize vector by its length
    Psi_d = atan(v[2],v[1])*180.0/pi #calculate sweep angle, convert to deg
    Theta_d = -asin(v[3])*180.0/pi #calculate theta angle, convert to deg

    return Psi_d,Theta_d
end

function createGeneralTransformationMatrix(angleArray,axisArray)
    #createGeneralTransformationMatrix  calculates general transformation matrix
    #   [dcmTotal] = createGeneralTransformationMatrix(angleArray,axisArray)
    #
    #   This function calculates the transformation matrix assocaited with a
    #   general Euler rotation sequence.
    #
    #      input:
    #      angleArray      = array of angles for Euler rotation sequence
    #      axisArray       = array of axis of rotatoins for Euler rotation
    #                        sequences
    #
    #      output:
    #      dcmTotal        = transformation matrix of specified euler rotation
    #                        sequence

    numRotations = length(angleArray) #get number of rotation to perform
    dcmArray = zeros(3,3,numRotations) #initialize individual rotation direction cosine matrix arrays

    for i=1:numRotations #calculate individual single rotatio direction cosine matrices
        dcmArray[:,:,i] = createSingleRotationDCM(angleArray[i],axisArray[i])
    end

    dcmTotal = dcmArray[:,:,1] #initialize dcmTotal as first rotation

    #multiply consecutive rotation sequency direction cosine matrices to arrive at overall transformation matrix
    for i=2:1:numRotations
        dcmTotal = dcmArray[:,:,i]*dcmTotal
    end

    return dcmTotal

end

function createSingleRotationDCM(angleDeg,axisNum)
    #This function creates a direction cosine matrix (dcm) associated
    #with a rotation of angleDeg about axisNum.

    angleRad = angleDeg*pi/180.0 #convert angle to radians

    if axisNum == 1 #direction cosine matrix about 1 axis
        dcm = [1.0 0.0 0.0
        0.0 cos(angleRad) sin(angleRad)
        0.0 -sin(angleRad) cos(angleRad)]
    elseif axisNum == 2 #direction cosine matrix about 2 axis
        dcm = [cos(angleRad) 0.0 -sin(angleRad)
        0.0 1.0 0.0
        sin(angleRad) 0.0 cos(angleRad)]
    elseif axisNum == 3 #direction cosine matrix about 3 axis
        dcm = [cos(angleRad) sin(angleRad) 0.0
        -sin(angleRad) cos(angleRad) 0.0
        0.0 0.0 1.0]
    else  #error catch
        error("Error: createSingleRotationDCM. Axis number must be 1, 2, or 3.")
    end

    return dcm

end

function calculateElementOrientation(mesh)
    #calculateElementOrientation calculates the orienation of elements in mesh
    #   [elOr] = calculateElementOrientation(mesh)
    #
    #   This function calculates the orientation of elements in a mesh.
    #
    #      input:
    #      mesh       = object containing mesh data
    #
    #      output:
    #      elOr       = object containing element orientation data

    numEl = mesh.numEl #get number of elements
    Psi_d=zeros(numEl) #initialize Psi, Theta, Twist, and Offset Arrays
    Theta_d=zeros(numEl)
    twist_d=zeros(numEl)
    Offset=zeros(3,numEl)    #offset is the hub frame coordinate of node 1 of the element
    elNum=zeros(numEl) #initialize element number array


    #calculate "mesh centroid"
    meshCentroid = [Statistics.mean(mesh.x) Statistics.mean(mesh.y) Statistics.mean(mesh.z)] #calculate a geometric centroid using all nodal coordinates
    lenv = zeros(numEl)
    for i = 1:numEl #loop over elements

        n1 = Int(mesh.conn[i,1]) #n1 := node number for node 1 of element i
        n2 = Int(mesh.conn[i,2]) #n2 := node number for node 2 of element i

        p1 = [mesh.x[n1] mesh.y[n1] mesh.z[n1]] #nodal coordinates of n1
        p2 = [mesh.x[n2] mesh.y[n2] mesh.z[n2]] #nodal coordinates of n2
        Offset[:,i] = p1 #set offset as position of n1

        v=p2-p1 #define vector from p1 to p2
        lenv[i] = LinearAlgebra.norm(v) #calculate element lengtt

        Psi_d[i],Theta_d[i] = calculatePsiTheta(v) #calculate elment Psi and Theta angles for orientation
        elNum[i] = mesh.conn[i,1] #get elemetn number

        nVec1,nVec2,nVec3 = rigidBodyRotation(0,0,1,[Psi_d[i],Theta_d[i]],[3,2]) #tranform a local normal "flapwise" vector in the element frame to the hub frame
        nVec = [nVec1 nVec2 nVec3]

        #for consistency, force the "flapwise" normal vector of an element to be
        #away from the machine

        # Mesh Type: 0-blade 1-tower 2-strut
        if mesh.type[i]==2
            refVector = [0;0;1]
        elseif mesh.type[i]==1
            refVector = [1;0;0]
        else
            refVector = p1-meshCentroid
        end

        refVector = refVector./LinearAlgebra.norm(refVector)
        dotTest = LinearAlgebra.dot(nVec,refVector)

        if dotTest<0 && abs(dotTest)>1.0e-4 #if vectors are not more or less aligned the normal vector is pointed inwards, towards the turbine (twist_d by 180 degrees)
            twist_d[i] = 180.0
        elseif abs(dotTest)<1.0e-4
            twist_dtemp = 90.0
            nVec1,nVec2,nVec3 = rigidBodyRotation(0,0,1,[Psi_d[i],Theta_d[i],twist_dtemp],[3,2,1])
            nVec = [nVec1 nVec2 nVec3]
            dotTest2 = LinearAlgebra.dot(nVec,refVector)
            if dotTest2<0 && abs(dotTest2)>1.0e-4 #if vectors are not more or less aligned the normal vector is pointed inwards, towards the turbine (twist_d by 180 degrees)
                twist_d[i] = twist_dtemp+180.0
                if abs(abs(twist_d[i])-270.0) < 1.0e-3
                    twist_d[i] = -90.0
                end
            else
                twist_d[i] = twist_dtemp
            end
        else  #the normal vector is pointed outwards, away from the turbine (no twist_d necessary)
            twist_d[i] = 0.0
        end

    end

    #assign data to element orientation (Ort) object
    return GyricFEA.Ort(Psi_d,Theta_d,twist_d,lenv,elNum,Offset)
end
