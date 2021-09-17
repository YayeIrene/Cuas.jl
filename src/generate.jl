

#function createCanon(u₀::Float64;hg=nothing, vg=nothing, zerotage=nothing,tc=nothing, lw=nothing,X2w=nothing, θₜ=0.0, ξₜ=0.0 )
 #    Canon(u₀,hg, vg, zerotage,tc, lw,X2w, θₜ, ξₜ )
#end
"""
    createTourelle(θ,ξ;optional arguments)

Creates an object turrets with elevation angle = θ and azimuth angle = ξ (in degrees).
optional arguments are angles (in degrees): Φ, Ψ, χ and rwY.
"""
function createTourelle(θ::Float64, ξ::Float64; Φ=nothing, ψ=nothing, χ=nothing, rWY=nothing)
    Tourelle(θ, ξ, Φ, ψ, χ, rWY)

end

#function createTarget(ρ::Float64)
#    Target(ρ)
#end
"""
    createProjectile(mass,Calibre;optional arguments)

Creates an object projectile with mass (kg), calibre (m).
optional arguments are: position, velocity, time of flight, moment of inertia, distance to the center of gravity.
"""
function createProjectile(mass::Float64, calibre::Float64; position=nothing, velocity =nothing, tof =nothing, Ix=nothing, Iy = nothing, Xcg=nothing)
    Projectile(mass, calibre, position, velocity, tof, Ix,Iy, Xcg)
end
