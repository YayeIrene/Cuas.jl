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
Returns the single shot hit probability (sshp) for a rectangular target
For the side of the rectangle standard deviations mean point of impact are defined in with "error" variable.
Vertical and horizontal errors are assumed to be independent...
```math
\\int_{-a}^{a}\\int_{-b}^{b}\\frac{1}{2\\pi\\sigma_x\\sigma_y}\\text{exp}[\\frac{(x-\\mu_x)^2}{2\\sigma_x^2}\\frac{(y-\\mu_y)^2}{2\\sigma_y^2}]dxdy
```
"""
function SSHP(target::TargetRect,error::RectErrorB)
    y = x->1/(sqrt(2*pi)*error.σ_x)*exp(-((x[1]-error.μ_x)^2)/(2*error.σ_x^2))*1/(sqrt(2*pi)*error.σ_y)*exp(-((x[2]-error.μ_y)^2)/(2*error.σ_y^2))
    g=hcubature(y,[-target.a,-target.b],[target.a, target.b])[1]
    return g
end

"""
    SSHP(target, error)

    Returns the single shot hit probability (sshp) for a circular target
    zero mean uncorrelated
    Vertical and horizontal errors are assumed to be independent...
```math
P(R)=\\int_{r=0}^{r=R} \\int_{\\theta=0}^{\\theta=2\\pi} \\frac{1}{2\\pi\\sigma_x\\sigma_y} \\text{exp}-[\\frac{(rcos\\theta)^2}{2\\sigma_x^2}+\\frac{(rsin\\theta)^2}{2\\sigma_y^2}]rdrd\\theta
```
"""
function SSHP(target::TargetCirc,error::RectErrorB)
    y=x->1/(2*pi*error.σ_x*error.σ_y)*exp(-((((x[1]*cos(x[2]))^2/(2*error.σ_x^2))+((x[1]*sin(x[2]))^2/(2*error.σ_y^2)))*x[1]))
    g=hcubature(y,[0,0],[target.ρ,2*pi])[1]
    return Pr
end

"""
    SSHP(target, error)

    Returns the single shot hit probability (sshp) for a spherical target
    zero mean uncorrelated
    Vertical and horizontal errors are assumed to be independent...
```math
P(R) = \\frac{1}{(2\\pi)^(3/2)\\sigma_x^2\\sigma_y^2\\sigma_z^2}\\int_{0}^{R}\\int_{0}^{\\pi}\\int_{0}^{2\\pi}f(r,\\phi,\\theta)r^2sin\\theta d\\theta d\\phi dr
```
    where
```
f(r,\\phi,\\theta) = exp-[\\frac{(r^2sin^2\\phi cos^2\\theta)}{2\\sigma_x^2}+\\frac{(r^2sin^2\\theta sin^2\\phi)}{2\\sigma_y^2}+\\frac{r^2cos^2\\phi}{2\\sigma_z^2}]rdrd\\theta
```
"""
function SSHP(target::TargetSphe,error::SpheError)
    y=x->1/((2*pi)^(3/2)*error.σ_x^2*error.σ_y^2*error.σ_z^2)*exp(-((((x[1]^2*(sin(x[3]))^2*(cos(x[2]))^2)/(2*error.σ_x^2))+((x[1]^2*(sin(x[2]))^2*(sin(x[3]))^2)/(2*error.σ_y^2))+(x[1]^2*(cos(x[3]))^2)/(2*error.σ_z^2))*x[1]^2*sin(x[2])))
    g=hcubature(y,[0,0,0],[target.ρ,2*pi,pi])[1]
    return Pr
end

"""
    MSHP(target,error,n)
Returns the multiple shot hit probability (mshp) of n rounds.
```math
P(n shots) = 1-(1-P)^n
```
"""
function MSHP(target::AbstractTarget,error::RectErrorB,n::Int64)
    Pi = SSHP(target,error)
    Pn = 1-(1-Pi)^n
    return Pn
end






"""
    CEP(σ)
Returns the circular error probable (CEP)
zero offset and equal variances (σ)
```math
CEP = 1.1774\\sigma
```
"""
function CEP(σ::Float64)
   CEP = 1.1774*σ
   return CEP
end

"""
    CEP(σx,σy)
Returns the circular error probable (CEP)
non equal variances (σx and σy)
```math
CEP \\cong \\sqrt{(\\sigma_x^2+\\sigma_y^2)\\big(1-\\frac{2(\\sigma_x^4+\\sigma_y^4)}{9(\\sigma_x^2+\\sigma_y^2)^2}\\big)^3}
```
"""
function CEP(σx::Float64, σy::Float64)
    CEP= sqrt((σx^2+σy^2)*(1-((2*(σx^4+σy^4)/(9*(σx^2+σy^2)^2))))^3)
    return CEP
end

"""
    CEP(error::RectErrorB;args)
Returns the circular error probable (CEP)
The mean point of impact can be different from (0,0).
Optinal arguments are:
* N the population of the sample size default value is 1e6
* p the probability default value is 50
"""
function CEP(error::RectErrorB;N=1000000, p=50)
mean = [error.μ_x,error.μ_y]
C = [error.σ_x^2 0; 0 error.σ_y^2]
d = MvNormal(mean, C)
x = rand(d, N)
s = Vector{Float64}(undef,N)
for i=1:N
    s[i] = sqrt(x[1,i]^2+x[2,i]^2)
end

s_o = sort(s)

CEP = s_o[Int(N*p/100)]
return CEP
end

"""
    SEP(σ)
Returns the circular error probable (CEP)
zero offset and equal variances (σ)
```math

```
"""
function SEP(σx::Float64, σy::Float64)

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
