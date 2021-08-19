#module types

#export Projectile

mutable struct RectErrorB
    σ_x::Float64
    σ_y::Float64
    μ_x::Float64
    μ_y::Float64
end

mutable struct Hull
    Φ::Union{Float64,Nothing}
    ψ::Union{Float64,Nothing}
    χ::Union{Float64,Nothing}
end

mutable struct Turret
     ξ::Union{Float64,Nothing}#azimuth
end

mutable struct Canon
     θ::Union{Float64,Nothing} #elevation
     lw::Float64
     u₀::Float64
end

 mutable struct Tourelle
     θ::Union{Float64,Nothing} #elevation
     ξ::Union{Float64,Nothing}#azimuth
     Φ::Union{Float64,Nothing}
     ψ::Union{Float64,Nothing}
     χ::Union{Float64,Nothing}
     rWY::Union{Float64,Nothing}
 end

 mutable struct Sight
     θ::Float64
     ξ::Float64
 end

mutable struct Tank
    hull::Hull
    turret::Turret
    canon::Canon
    sight::Sight
    rWY::Float64
end

 #mutable struct Target
#     ρ::Float64
# end
