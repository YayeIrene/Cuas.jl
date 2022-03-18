#module types

#export Projectile

mutable struct RectErrorB <:AbstractTarget
    σ_x::Float64
    σ_y::Float64
    μ_x::Float64
    μ_y::Float64
end



 mutable struct Tourelle
     θ::Union{Float64,Nothing} #elevation
     ξ::Union{Float64,Nothing}#azimuth
     Φ::Union{Float64,Nothing}
     ψ::Union{Float64,Nothing}
     χ::Union{Float64,Nothing}
     rWY::Union{Float64,Nothing}
 end

 """
 Target(θ,ϕ)
defines a target
* t0 : the time at which the target appears
* ρ : the range to the target in m
* θ : the elevation of the LOS
* ϕ : the azimuth of the LOS
* ρ′ : target velocity in range
* θ′ : target velocity in elevation
* ϕ′ : target velocity in azimuth
* ρ″ : target acceleration in range
* θ″ : target acceleration in elelvation
* ϕ″ : target accelertion in azimuth
* type : type of target
"""
 mutable struct Target <:AbstractTarget
    t0::Float64
    ρ::Float64
    θ::Float64
    ϕ::Float64
    ρ′::Float64
    θ′::Float64
    ϕ′::Float64
    ρ″::Float64
    θ″::Float64
    ϕ″::Float64
    type::String
    size::Union{Vector{Float64},Nothing}
    alive::Bool
  function Target(θ::Float64, ϕ::Float64)
      new(0.0, 0.0,θ,ϕ,0.0,0.0,0.0,0.0,0.0,0.0,"unkown",nothing,true)
  end
end

 mutable struct Uav <:AbstractTarget 
    ρ::Union{Float64,Nothing} 
    position::Union{Vector{Float64},Nothing}
    velocity::Union{Vector{Float64},Nothing}
    θ::Union{Float64,Nothing}
    ϕ::Union{Float64,Nothing}
    alive::Bool
    function Uav(position::Vector{Float64},velocity::Vector{Float64})
        new(nothing,position,velocity,nothing,nothing,true)
    end
 end 

 #mutable struct Target
#     ρ::Float64
# end
