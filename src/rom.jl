
"""

  reducedOrderModel(elStorage,feamodel,mesh,el,displ)

This function executes a reduced order model analysis.

#Input
* `elStorage`      object containing stored element matrices
* `feamodel`          object containing feamodel information
* `mesh`           object containing mesh information
* `el`             object containing element information
* `displ`          displacement vector for use in pre-stressed analysis

#Output
* `rom`            object containing a reduced order feamodel
"""
function reducedOrderModel(elStorage,feamodel,mesh,el,displ)

	countedNodes = []

	rom0 = calculateROM(feamodel,mesh,el,displ,zeros(9),elStorage,countedNodes)  #Omega   = 0 #calculates system matrices for parked condition

	#calculates system matrices for various acceleration, rotor speed, rotor accelration combinations
	omx = 1
	omy = 1
	omz = 1

	omxdot = 1
	omydot = 1
	omzdot = 1

	a_x = 1
	a_y = 1
	a_z = 1

	rom1 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 omx 0.0 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_x = 1 #OmegaDot_i = 0 accel _i = 0
	rom2 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 omy 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_y = 1 #OmegaDot_i = 0 accel_i = 0
	rom3 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 0.0 omz 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_z = 1 #OmegaDot_i = 0 accel_i = 0
	rom4 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 omx omy 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_x,y = 1 #OmegaDot_i = 0 accel_i = 0
	rom5 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 omy omz 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_y,z = 1 #OmegaDot_i = 0 accel_i = 0
	rom6 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 omx 0.0 omz 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_x,z = 1 #OmegaDot_i = 0 accel_i = 0
	rom7 = calculateROMGyric(feamodel,mesh,el,displ,[a_x 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_i = 0 accel_x = 1
	rom8 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 a_y 0.0 0.0 0.0 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_i = 0 accel_y = 1
	rom9 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 a_z 0.0 0.0 0.0 0.0 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_i = 0 accel_z = 1
	rom10 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 0.0 0.0 omxdot 0.0 0.0],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_1 = 0 accel_i=0
	rom11 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 0.0 0.0 0.0 omydot 0.0],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_2 = 0 accel_i=0
	rom12 = calculateROMGyric(feamodel,mesh,el,displ,[0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 omzdot],elStorage,rom0,countedNodes)  #Omega_i = 0 #OmegaDot_3 = 0 accel_i=0

	#reduced order structural stiffness, mass, and damping
	Kr = rom0.Kr
	Mr = rom0.Mr
	Cr = rom0.Cr + rom0.CrModal
	Fr = rom0.Fr

	#reduced order transformation matrices
	Phi = rom0.Phi
	invPhi = rom0.invPhi

	#reduced order spin softening coefficient matrices
	SrOx2 = rom1.Kr
	SrOy2 = rom2.Kr
	SrOz2 = rom3.Kr
	SrOxOy = rom4.Kr - rom1.Kr - rom2.Kr
	SrOyOz = rom5.Kr - rom2.Kr - rom3.Kr
	SrOxOz = rom6.Kr - rom1.Kr - rom3.Kr

	#reduced order centrifugal load coefficient matrices
	FrOx2 = rom1.Fr
	FrOy2 = rom2.Fr
	FrOz2 = rom3.Fr
	FrOxOy = rom4.Fr - rom1.Fr - rom2.Fr
	FrOyOz = rom5.Fr - rom2.Fr - rom3.Fr
	FrOxOz = rom6.Fr - rom1.Fr - rom3.Fr
	FrOxdot = rom10.Fr
	FrOydot = rom11.Fr
	FrOzdot = rom12.Fr

	#need to calculate reduced order acceleration/body force coefeficients
	FrAx = rom7.Fr
	FrAy = rom8.Fr
	FrAz = rom9.Fr

	#reduced order gyric coefficient matrices
	GrOx = rom1.Cr
	GrOy = rom2.Cr
	GrOz = rom3.Cr

	#reduced order circulatory coefficient matrices
	HrOx = 0.5*rom1.Cr
	HrOy = 0.5*rom2.Cr
	HrOz = 0.5*rom3.Cr

	rom = ROM(Kr,Mr,Cr,0.0,Fr,Phi,invPhi,SrOx2,SrOy2,SrOz2,SrOxOy,SrOyOz,SrOxOz,FrOx2,FrOy2,FrOz2,FrOxOy,FrOyOz,FrOxOz,FrOxdot,FrOydot,FrOzdot,FrAx,FrAy,FrAz,GrOx,GrOy,GrOz,HrOx,HrOy,HrOz)

	return rom
end

"""

	calculateROM(model,mesh,el,displ,omegaVec,omegaDotVec,elStorage,countedNodes)

This function calculates a reduced order model for a conventional
structural dynamics system (parked, non-rotating)

#Input
* `model`        object containing model data
* `mesh`         object containing mesh data
* `el`           object containing elementdata
* `displ`        displacement vector
* `rbData`:     vector containing rigid body displacement, velocity, and acceleration
* `elStorage`    object containing stored element data
* `countedNodes`    prevents applied nodal terms from double counting

#Output
* `rom`          object containing reduced order model data
"""
function  calculateROM(feamodel,mesh,el,displ,rbData,elStorage,countedNodes)
	feamodel.analysisType = "ROM"
	BC = feamodel.BC

	elementOrder = feamodel.elementOrder
	numNodesPerEl = elementOrder + 1
	numDOFPerNode = 6
	totalNumDOF = mesh.numNodes * numDOFPerNode

	eldisp = zeros(numNodesPerEl*numDOFPerNode)
	Kg = zeros(totalNumDOF,totalNumDOF)
	Mg = zeros(totalNumDOF,totalNumDOF)
	Cg = zeros(totalNumDOF,totalNumDOF)
	Fg = zeros(totalNumDOF)

	TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,displ,0.0,elStorage; Kg,Mg,Cg,Fg,rbData,countedNodes)

	#----------------------------------------------------------------------
	#APPLY CONSTRAINT
	Kg = applyConstraints(Kg,feamodel.jointTransform)   #apply constraints to system matrices and force vector
	Mg = applyConstraints(Mg,feamodel.jointTransform)
	Cg = applyConstraints(Cg,feamodel.jointTransform)
	Fg = applyConstraintsVec(Fg,feamodel.jointTransform)

	#APPLY BOUNDARY CONDITIONS
	KgTotal = applyBCModal(Kg,feamodel.BC.numpBC,BC.map) #apply boundary conditions to system matrices and force vector
	MgTotal = applyBCModal(Mg,feamodel.BC.numpBC,BC.map)
	CgTotal = applyBCModal(Cg,feamodel.BC.numpBC,BC.map)
	FgTotal = applyBCModalVec(Fg,feamodel.BC.numpBC,BC.map)

	sysMat=MgTotal\KgTotal                      #eigenvalue solve on spring mass system only
	nev = size(sysMat)[1]#min(size(sysMat)[1],feamodel.numModes)
    eigVal, eigVec = ArnoldiMethod.partialeigen(ArnoldiMethod.partialschur(sysMat; nev, which=ArnoldiMethod.LM())[1])

	perm = sortperm(eigVal, by=(eigVal)->(abs(eigVal),imag(eigVal)), rev=false)
    # eigVal .= eigVal[perm] #NOT USED
    eigVec .= eigVec[:,perm]
	eigVec = eigVec[:,1:feamodel.numModes] # shorten to the required length

	# for icol = 1:feamodel.numModes
	# 	if eigVec[1,icol] < 0.0
	# 		eigVec[1,icol] .*= -1
	# 	end
	# end

	len = length(KgTotal)
	Phi = real(eigVec)   #assign modal matrix to Phi
	invPhi = LinearAlgebra.pinv(Phi) #calculate psuedo inverse of Phi

	Mr = Phi'*MgTotal*Phi   #calculate ROM mass matrix
	Kr = Phi'*KgTotal*Phi   #calculate ROM stiffness matrix
	Cr = Phi'*CgTotal*Phi   #calculate ROM damping matrix
	Fr = Phi'*FgTotal       #calculate ROM force vector

	#apply modal damping
	dampRatio = 0.02
	Cmodal = zero(Kr)
	for i=1:length(Kr[:,1])
		Cmodal[i,i] = 2*dampRatio*sqrt(Kr[i,i]*Mr[i,i])
	end
	CrModal = Cmodal

	return ROM(Kr,Mr,Cr,CrModal,Fr,Phi,invPhi)
end


#TODO: update applyBCModal to detect if the K is a vector or matrix and do this (since it is 99# the same)
function applyBCModalVec(F,numpBC,bcMap)
	#applyBCModal Applies boundary conditions to system for modal analysis
	#   [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
	#
	#   This function applies boundary conditions to a system matrix for modal
	#   analysis
	#
	#      #Input
	#      K             = assembled global system matrix
	#      BC            = struct of boundary condition information
	#      numDofPerNode = number of degrees of freedom per node

	#      #Output
	#      K             = global system matrix with boundary conditions
	#      dofVector     = reduced DOF vector after imposing BCs

	numEq=size(F)[1]
	Fnew = zeros(numEq - numpBC)
	indVec = zeros(Int,numEq - numpBC)

	index = 1
	for i=1:numEq
		if bcMap[i] != -1
			indVec[index] = i
			index = index +1
		end
	end

	#APPLY BCs FOR PRIMARY VARIABLE
	if numpBC > 0
		Fnew = F[indVec]
	else
		Fnew = F
	end
	return Fnew
end

"""

	structuralDynamicsTransientROM(feamodel,mesh,el,dispData,Omega,OmegaDot,time,delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData)

Performs transient structural dynamics analysis using a reduced order feamodel (ROM).

#Input
* `feamodel`:      object containing feamodel data
* `mesh`:       object containing mesh data
* `el`:         object containing element data
* `dispData`:   object containing displacement data
* `Omega`:      rotor speed (Hz)
* `OmegaDot`:   rotor acceleratin (Hz)
* `time`:       current simulation time
* `delta_t`:    time step size
* `elStorage`:  object containing stored element data
* `rom`:        object containing reduced order feamodel represnetation
* `Fexternal`:  vector containing external force values
* `Fdof`:       vector containing global DOF numbering associated with external force values
* `CN2H`:       transformation matrix from inertial frame to hub frame
* `rbData`:     vector containing rigid body displacement, velocity, and acceleration

#Output
* `dispOut`:       object containing displacement data at end of time step
* `FReaction_sp1`: vector containing reaction force at turbine base at end of time step
"""
function  structuralDynamicsTransientROM(feamodel,mesh,el,dispData,Omega,OmegaDot,time,delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData)
	feamodel.analysisType = "ROM" #enforce analysis type based on function call
	###-------- get feamodel information -----------
	numEl = mesh.numEl
	x = mesh.x
	y = mesh.y
	z = mesh.z
	conn = Int.(mesh.conn)
	numNodes = length(x)
	elementOrder = feamodel.elementOrder
	BC = feamodel.BC
	nlROM = feamodel.nlOn

	numNodesPerEl = elementOrder + 1
	numDOFPerNode = 6
	totalNumDOF = numNodes * numDOFPerNode
	_,numReducedDOF = size(feamodel.jointTransform)
	nodalTerms = feamodel.nodalTerms
	countedNodes = []
	###-----------------------------------------

	###------- intitialization -----------------
	###-------------------------------------------
	analysisType = "M"
	###------ newmark integration parameters ---------
	alpha = 0.5
	gamma = 0.5
	beta = 0.5*gamma

	delta_t = delta_t
	a1 = alpha*delta_t
	a2 = (1.0-alpha)*delta_t
	a3 = 1.0/(beta*delta_t*delta_t)
	a4 = a3*delta_t
	a5 = 1.0/gamma-1.0
	a6 = alpha/(beta*delta_t)
	a7 = alpha/beta - 1.0
	a8 = delta_t*(alpha/gamma-1.0)
	timeInt = TimeInt(delta_t,a1,a2,a3,a4,a5,a6,a7,a8)
	disp_s = dispData.displ_s
	dispdot_s = dispData.displdot_s
	dispddot_s = dispData.displddot_s
	###-----------------------------------------------

	#initialize displacements, tolerance, uNorm, iteration count for nonlinear
	#iteration
	displ_iter = copy(disp_s)
	displ_last = copy(disp_s)
	iterationCount = 1
	maxIter = 50
	uNorm = 1.0e6
	tol = 1.0e-6

	elx = zeros(numNodesPerEl) #initialize element coordinate list

	eta_iter = zeros(feamodel.numModes)
	etadot_iter = zeros(feamodel.numModes)
	etaddot_iter = zeros(feamodel.numModes)

	while iterationCount < maxIter && uNorm > tol  #iteration loop
		Kg = zeros(totalNumDOF,totalNumDOF)  #initialize global stiffness and force vector
		Fg = zeros(totalNumDOF)
		## Element Assembly Loop for NL Terms
		if nlROM
			for i=1:numEl
				eldispiter = zeros(numNodesPerEl*numDOFPerNode)
				#Calculate Ke and Fe for element i
				index = 1                           #initialize element data
				analysisType = analysisType
				elementOrder = elementOrder
				modalFlag = true
				xloc = [0.0 el.elLen[i]]
				sectionProps = el.props[i]
				sweepAngle = el.psi[i]
				coneAngle = el.theta[i]
				rollAngle = el.roll[i]
				aeroSweepAngle = 0.0

				for j=1:numNodesPerEl

					elx[j] = x[conn[i,j]]

					#get element nodal displacements at s and s-1 time step
					for k=1:numDOFPerNode
						eldispiter[index] = displ_iter[(conn[i,j]-1)*numDOFPerNode + k]
						index = index + 1
					end
				end

				disp= eldispiter
				omegaVec = zeros(3) #convert from platform frame to hub-frame
				omegaDotVec = zeros(3)
				Omega = 0.0
				OmegaDot = 0.0
				CN2H = CN2H
				preStress = false
				useDisp = true
				iterationType = "DI"  #uses direct iteration

				#do element calculations
				elInput = ElInput(elementOrder,modalFlag,timeInt,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,iterationType,useDisp,preStress,false,false,1.0,1.0,1.0,1.0,1.0,analysisType,disp,zero(disp),zero(disp),zero(disp),0.0,0.0,0.0,0.0,elx,0.0,0.0,true,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,false)
				Ke = calculateTimoshenkoElementNLSS(elInput) #calculate nonlinear timoshenko element stiffness matrix

				Kg = assemblyMatrixOnly(Ke,conn[i,:],numNodesPerEl,numDOFPerNode,Kg) #assemble nonlinear timoshenko element stiffness matrix
				#................................................
			end
		end
		###------- end element calculation and assembly ------------------
		##
		
		#Apply external loads to structure
		for i=1:length(Fexternal)
			Fg[Fdof[i]] = Fg[Fdof[i]] + Fexternal[i]
		end

		###------ apply constraints on system -----------------------------------
		##
		Fg = applyConstraintsVec(Fg,feamodel.jointTransform)

		if nlROM
			Kg = applyConstraints(Kg,feamodel.jointTransform)
		end
		###----------------------------------------------------------------------
		##
		#Apply BCs to global system
		Fg     = applyBCModalVec(Fg,BC.numpBC,BC.map)

		eta_s = dispData.eta_s
		etadot_s = dispData.etadot_s
		etaddot_s = dispData.etaddot_s

		Phi       = rom.Phi
		Feta = (Phi')*Fg   #transform global displacement vector to modal space

		if nlROM
			Kgnl   = applyBCModal(Kg,BC.numpBC,BC.map)
			Kgnlrom = Phi'*Kgnl*Phi     #transform nonlinear stiffness matrix to modal space
		else
			Kgnlrom = zeros(size(rom.Mr)) #if nonlinear deactivated set Kgnlrom to zero matrix
		end


		#define omega_i and omegaDot_i and body accelerations
		omega_platform_s = rbData[4:6]
		omega_x = omega_platform_s[1]
		omega_y = omega_platform_s[2]
		omega_z = omega_platform_s[3] + Omega*2*pi

		omega_platform_dot = rbData[7:9]
		omegaDot_x = omega_platform_dot[1]
		omegaDot_y = omega_platform_dot[2]
		omegaDot_z = omega_platform_dot[3] + OmegaDot*2*pi

		a_x = rbData[1] #platform accelerations (in hub frame)
		a_y = rbData[2]
		a_z = rbData[3]

		if feamodel.gravityOn
			g= 9.81 #gravitational acceleration m/s^2
		else
			g = 0.0
		end

		a_x_n = 0.0 #accelerations in inertial frame
		a_y_n = 0.0
		a_z_n = g
		a_temp = CN2H*[a_x_n; a_y_n; a_z_n]

		a_x = a_x + a_temp[1]
		a_y = a_y + a_temp[2]
		a_z = a_z + a_temp[3]

		#calculate reduced order spin soft, cent force, body force,
		#circulatory, coriolis

		Seff = rom.SrOx2.*omega_x^2 + rom.SrOy2.*omega_y^2 + rom.SrOz2.*omega_z^2 +
		rom.SrOxOy.*omega_x*omega_y + rom.SrOyOz.*omega_y*omega_z +
		rom.SrOxOz.*omega_x*omega_z

		FcentEff = rom.FrOx2.*omega_x^2 + rom.FrOy2.*omega_y^2 + rom.FrOz2.*omega_z^2 +
		rom.FrOxOy.*omega_x*omega_y + rom.FrOyOz.*omega_y*omega_z +
		rom.FrOxOz.*omega_x*omega_z + rom.FrOxdot*omegaDot_x +
		rom.FrOydot*omegaDot_y + rom.FrOzdot*omegaDot_z

		FbodyEff = rom.FrAx.*a_x + rom.FrAy.*a_y + rom.FrAz.*a_z

		Heff = rom.HrOx.*omegaDot_x + rom.HrOy.*omegaDot_y + rom.HrOz.*omegaDot_z

		Geff = rom.GrOx*omega_x + rom.GrOy*omega_y + rom.GrOz*omega_z

		#combine for effective reduced order stiffness, damping, force
		Keff = rom.Kr + Seff + Heff + Kgnlrom
		Ceff = rom.Cr + Geff
		Feff = Feta + FcentEff + FbodyEff

		eta_iter,etadot_iter,etaddot_iter = timeIntegrateSubSystemEff(rom.Mr,Keff,Ceff,Feff,timeInt,eta_s,etadot_s,etaddot_s)

		#reconstruct constrained dof vector from boundary conditions
		dispVec = Phi*eta_iter
		dispdotVec = Phi*etadot_iter
		dispddotVec = Phi*etaddot_iter

		#[dispVec,dispdotVec,dispddotVec]    = constructReducedDispVector(dispVec,dispdotVec,dispddotVec,numNodes,BC)
		dispVec = constructReducedDispVectorSingle(dispVec,BC.redVectorMap)

		displ_iter = feamodel.jointTransform*dispVec #TODO: mutate not overwrite

		if nlROM
			uNorm = calcUnorm(displ_iter,displ_last)
			displ_last = copy(displ_iter)
			iterationCount = iterationCount + 1
		else
			uNorm = 0.0
		end

	end
	if iterationCount>=maxIter
		error("Maximum iterations exceeded.")
	else

	end
	###------ calculate reaction at turbine base ----------------------------
	reactionNodeNumber = feamodel.platformTurbineConnectionNodeNumber
	FReaction_sp1 = calculateReactionForceAtNode(reactionNodeNumber,feamodel,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H,countedNodes)
	###----------------------------------------------------------------------
	#Calculate strain
	elStrain = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ_iter,feamodel.nlOn)

	dispdotVec    = constructReducedDispVectorSingle(dispdotVec,BC.redVectorMap)  # reconstruct reduced velocity vector with boundary conditions
	dispddotVec    = constructReducedDispVectorSingle(dispddotVec,BC.redVectorMap) # reconstruct reduced acceleration vector with boundary conditions
	displdot_iter = feamodel.jointTransform*dispdotVec   #construct unconstrained velocity vector
	displddot_iter = feamodel.jointTransform*dispddotVec  #construct unconstrained acceleration vector

	dispOut = DispOut(elStrain,displ_iter,displddot_iter,displdot_iter,eta_iter,etadot_iter,etaddot_iter)

	return elStrain,dispOut,FReaction_sp1

end

function constructReducedDispVectorSingle(vec1,redVectorMap)
	#This function reconstructs displacement vector to account for boundary
	#conditions
	len = length(redVectorMap)
	vec1Red = zeros(len)

	for i=1:len
		if redVectorMap[i] == -1.0
			vec1Red[i] = 0.0
		else
			vec1Red[i] = vec1[Int(redVectorMap[i])]
		end
	end
	return vec1Red
end

"""

  calculateROMGyric(feamodel,mesh,el,displ,omegaVec,omegaDotVec,elStorage,rom0,countedNodes)

Calculates a reduced order feamodel with rotational/ rigid
body motion effects

#Input
* `feamodel`:        object containing feamodel data
* `mesh`:         object containing mesh data
* `el`:           object containing elementdata
* `displ`:        displacement vector
* `rbData`:     vector of hub frame accel (1-3), angular velocity components (4-6), and angular accleration (7-9)
* `elStorage`:    object containing stored element data
* `rom0`:         object containing parked/conventional reduced order feamodel
* `countedNodes`:    prevents applied nodal terms from double counting

#Output
* `rom`:          object containing reduced order feamodel data
"""
function  calculateROMGyric(feamodel,mesh,el,displ,rbData,elStorage,rom0,countedNodes)
	feamodel.analysisType = "ROM" #force analysis type to match function call

	elementOrder = feamodel.elementOrder
	numNodesPerEl = elementOrder + 1
	numDOFPerNode = 6
	totalNumDOF = mesh.numNodes * numDOFPerNode
	Kg = zeros(totalNumDOF,totalNumDOF)
	Cg = zeros(totalNumDOF,totalNumDOF)
	Fg = zeros(totalNumDOF)
	eldisp = zeros(numNodesPerEl*numDOFPerNode)

	TimoshenkoMatrixWrap!(feamodel,mesh,el,eldisp,displ,0.0,elStorage;
	    Kg,Cg,Fg,rbData,countedNodes)

	###----------------------------------------------------------------------
	#APPLY CONSTRAINT
	Kg = applyConstraints(Kg,feamodel.jointTransform) #apply constraints to system matrices and force vector
	Cg = applyConstraints(Cg,feamodel.jointTransform)
	Fg = applyConstraintsVec(Fg,feamodel.jointTransform)


	#APPLY BOUNDARY CONDITIONS
	KgTotal = applyBCModal(Kg,feamodel.BC.numpBC,feamodel.BC.map) #apply boundary conditions to system matrices and force vector
	CgTotal = applyBCModal(Cg,feamodel.BC.numpBC,feamodel.BC.map)
	FgTotal = applyBCModalVec(Fg,feamodel.BC.numpBC,feamodel.BC.map)

	Phi = rom0.Phi

	Kr = Phi'*KgTotal*Phi - rom0.Kr #removes contribution of structural stiffness to reduced K
	Cr = Phi'*CgTotal*Phi - rom0.Cr #removes contribution of structural damping to reduced C
	Fr = Phi'*FgTotal - rom0.Fr

	return ROM(Kr,Cr,Fr)

end

"""

  timeIntegrateSubSystemEff(M,K,C,F,timeInt,u,udot,uddot)

Performs integration of a system using the Newmark-Beta
method(constant-average acceleration sceheme). The integration
parameters are calculated before hand and store in the timeInt object.

#Input
* `M`        system mass matrix
* `K`        system sttiffness matrix
* `C`        system damping matrix
* `F`        system force vector
* `timeInt`  object containing time integraton parameters
* `u`        displacement at beginning of time step
* `udot`     velocity at beginning of time step
* `uddot`    acceleration at beginning of time step


#Output
* `unp1`:        displacement at end of time step
* `udotnp1`:     velocity at end of time step
* `uddotnp1`:    acceleration at end of time step
"""
function timeIntegrateSubSystemEff(M,K,C,F,timeInt,u,udot,uddot)

	#transient integration of sub system using Newmark Beta method
	a1 = timeInt.a1
	a2 = timeInt.a2
	a3 = timeInt.a3
	a4 = timeInt.a4
	a5 = timeInt.a5
	a6 = timeInt.a6
	a7 = timeInt.a7
	a8 = timeInt.a8
	A = a3*u + a4*udot + a5*uddot
	B = a6*u + a7*udot + a8*uddot

	Khat = K + a3.*M + a6.*C
	Fhat = F + M*(A) + C*(B)

	unp1 = Khat\Fhat

	uddotnp1 = a3*(unp1-u) - a4*udot - a5*uddot
	udotnp1 =  udot + a2*uddot + a1*uddotnp1

	return unp1,udotnp1,uddotnp1

end
