module Cuas
using ExternalBallistics
using HCubature
using ReferenceFrameRotations
# Write your package code here.
include("types.jl")
include("roundPerRound.jl")
include("generate.jl")
include("rotations.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, Canon,createTarget, createTourelle, createCanon, createProjectile, targetPos, muzzlePos, muzzleVel, wind


end
