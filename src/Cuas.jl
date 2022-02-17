module Cuas
using ExternalBallistics
using WeaponSystems
using ErrorBudget
using HCubature
using ReferenceFrameRotations, Distances
using Distributions, DataFrames
using Vulnerability, OpticSim
using SimJulia, ResumableFunctions

# Write your package code here.
include("types.jl")
include("roundPerRound.jl")
include("generate.jl")
include("rotations.jl")
include("aimpoint.jl")
include("burst.jl")
include("sensitivity.jl")
include("firingChain.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, createTourelle, createProjectile,  Tourelle,
 SSHP,CEP, Target, phitAbm,phitAbmBurst, impactAbm, phit, sshpAbm, acquisition


end
