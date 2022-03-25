module Cuas
using ExternalBallistics
using WeaponSystems
using ErrorBudget
using HCubature
using ReferenceFrameRotations, Distances,LinearAlgebra
using Distributions, DataFrames
using Vulnerability, OpticSim
using SimJulia, ResumableFunctions
using CoordinateTransformations

# Write your package code here.
include("types.jl")

include("generate.jl")
include("rotations.jl")
include("aimpoint.jl")
include("burst.jl")
include("sensitivity.jl")
include("ballisticComputer.jl")
include("firingChain.jl")
include("abm.jl")
include("roundPerRound.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, createTourelle, createProjectile,  Tourelle,
 SSHP,CEP, Target, phitAbm,phitAbmBurst, impactAbm, phit, sshpAbm, acquisition, hit,
 phitAtLeastOneBurst, Uav, aiming!


end
