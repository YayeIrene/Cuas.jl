module Cuas
using ExternalBallistics
using WeaponSystems
using ErrorBudget
using HCubature
using ReferenceFrameRotations
using Distributions
# Write your package code here.
include("types.jl")
include("roundPerRound.jl")
include("generate.jl")
include("rotations.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, createTourelle, createProjectile,  Tourelle,
 SSHP,CEP


end
