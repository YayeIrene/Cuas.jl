module Cuas
using ExternalBallistics
using HCubature
# Write your package code here.
include("types.jl")
include("roundPerRound.jl")
#include("phit.jl")
#using .phit

export pHit, RectErrorB, Gun


end
