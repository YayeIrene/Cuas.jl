

function createCanon(u₀::Float64;hg=nothing, vg=nothing, zerotage=nothing,tc=nothing, lw=nothing,X2w=nothing, θₜ=0.0, ξₜ=0.0 )
    Canon(u₀,hg, vg, zerotage,tc, lw,X2w, θₜ, ξₜ )
end
function createTourelle(θ::Float64, ξ::Float64; Φ=nothing, ψ=nothing, χ=nothing, rWY=nothing)
    Tourelle(θ, ξ, Φ, ψ, χ, rWY)

end

function createTarget(ρ::Float64)
    Target(ρ)
end

function createProjectile(mass::Float64, calibre::Float64; position=nothing, velocity =nothing, tof =nothing, Ix=nothing, Iy = nothing, Xcg=nothing)
    Projectile(mass, calibre, position, velocity, tof, Ix,Iy, Xcg)
end
