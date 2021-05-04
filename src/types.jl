module types

export Projectile

mutable struct Projectile
    position::Float64
    velocity :: Float64
    mass :: Float64
    calibre :: Float64
end

mutable struct Fragment
    position::Float64
    velocity :: Float64
    mass :: Float64
    calibre :: Float64
end

mutable struct Target #momenteel niet gebruikt
    θ :: Float64 # angle w.r.t. artificial north
    ρ :: Float64 # distance to own ship
    v :: Float64 # speed of the target
    d :: Float64 # diameter of the target
    A :: Float64 #target area
    Type :: Int64 #0: UAV, 1: Missile, 2:FIAC/LFA

end

mutable struct Results
    veject::Array{Float64,1}
    pvelocity::Array{Float64,1}
    coneangle::Array{Float64,1}
    alethal::Array{Float64,1}
    rlethal::Array{Float64,1}
    nbrfragment::Array{Float64,1}
    vifragment::Array{Float64,1}
    vffragment::Array{Float64,1}
    sshp::Array{Float64,1}
    pkfragment::Array{Float64,1}
    pkprojectile::Array{Float64,1}
    pkh::Array{Float64,1}
    rounds::Float64
    missed::Float64
    ekin::Float64
    pk::Float64
end

end
