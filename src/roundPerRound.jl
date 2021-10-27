#module roundPerRound

function σtot(σ::Float64,ρ::Float64)
    ρ*σ*1e-3
end

function SSHP(projectile::AbstractPenetrator, target::TargetCirc)
    D = abs(projectile.position-target.position)
    Stot=σtot(projectile.σ,D)
    y = (x)->x[1]/Stot^2*exp(-0.5*x[1]^2/Stot^2)
    g=hcubature(y,0,0.5*target.ρ)[1]
    return g
end
"""
    SSHP(target,error)
Returns the single shot hit probability (sshp) for a rectangular target.
"""
function SSHP(target::TargetRect,error::RectErrorB)
    y = x->1/(sqrt(2*pi)*error.σ_x)*exp(-((x[1]-error.μ_x)^2)/(2*error.σ_x^2))*1/(sqrt(2*pi)*error.σ_y)*exp(-((x[2]-error.μ_y)^2)/(2*error.σ_y^2))
    g=hcubature(y,[-target.a,-target.b],[target.a, target.b])[1]
    return g
end



"""
    pHit(target,error)

Returns the probability of hitting the target defined by'target' if the estimated error is provided by 'error'.
```math
\\int_{-a}^{a}\\int_{-b}^{b}\\frac{1}{2\\pi\\sigma_x\\sigma_y}\\text{exp}[\\frac{(x-\\mu_x)^2}{2\\sigma_x^2}\\frac{(y-\\mu_y)^2}{2\\sigma_y^2}]dxdy
```
"""
function pHit(target::AbstractTarget, error::RectErrorB)
    SSHP(target, error)
end
#end
#\\[\\frac{(x-\\mu_x)^2}{2\\sigma_x^2}
#+\\frac{(y-\\mu_y)^2}{2\\sigma_y^2}\\]dxdy\$
