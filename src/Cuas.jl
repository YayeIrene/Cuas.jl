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

export pHit, RectErrorB, createTourelle, createProjectile, targetPos, muzzlePos, muzzleVel, wind, Tourelle,
Hull, Turret, Canon, Sight, Tank, SSHP


end
