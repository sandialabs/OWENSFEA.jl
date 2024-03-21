module OWENSFEA

import LinearAlgebra
import SparseArrays
import ArnoldiMethod
import DelimitedFiles
import Printf
import MAT
import StaticArrays

export mapACloads, calculateStructureMassProps, createJointTransform, calculateReducedDOFVector, constructReducedDispVectorMap, calculateBCMap, calculateElementMass!, ConcMassAssociatedWithElement, applyBC, calculateReactionForceAtNode, calculateStrainForElements, findElementsAssociatedWithNodeNumber, applyGeneralConcentratedTerms


include("structs.jl")
include("unsteady.jl")
include("modal.jl")
include("utilities.jl")
include("timoshenko.jl")
include("steady.jl")
include("rom.jl")
include("intermediate.jl")

end
