#module types

#export Projectile

mutable struct RectErrorB
    σ_x::Float64
    σ_y::Float64
    μ_x::Float64
    μ_y::Float64
end

mutable struct Gun
    hg::Float64
    vg::Float64
    zerotage::Float64
end
