module roundPerRound
using ExternalBallistics
using QuadGK

export pHit

function σtot(σ::Float64,ρ::Float64)
    ρ*σ*1e-3
end

function SSHP(projectile::Projectile1D, target::Target1D)
    D = abs(projectile.position-target.position)
    Stot=σtot(projectile.σ,D)
    y = (x)->x/Stot^2*exp(-0.5*x^2/Stot^2)
    g=quadgk(y,0,0.5*target.size)[1]
    return g
end


function pHit(projectile::Projectile1D,target::Target1D;HIT_PROBABILITY=0.95, MAX_BURST_SIZE = 150)
    Pnhit=1.0
    burst=0

while (1-Pnhit) < HIT_PROBABILITY && burst<MAX_BURST_SIZE
        Phiti =  SSHP(projectile,target)
        Pnhit=Pnhit*(1-Phiti)
        burst=burst+1

end
return 1-Pnhit, burst

end
end
