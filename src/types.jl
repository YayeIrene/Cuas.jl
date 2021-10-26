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

 

 #mutable struct Target
#     ρ::Float64
# end
