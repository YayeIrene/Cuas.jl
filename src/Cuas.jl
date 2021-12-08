module Cuas
using ExternalBallistics
using WeaponSystems
using ErrorBudget
using HCubature
using ReferenceFrameRotations, Distances
using Distributions, DataFrames
using Vulnerability, OpticSim
# Write your package code here.
include("types.jl")
include("roundPerRound.jl")
include("generate.jl")
include("rotations.jl")
include("aimpoint.jl")
include("burst.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, createTourelle, createProjectile,  Tourelle,
 SSHP,CEP, Target,adjustedFire!, phitAbm,phitAbmBurst


end
