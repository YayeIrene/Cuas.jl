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

"""
phitAbm(tank,t, proj, w,atmosphere, fuze, aero,zones, shapes,dfError,linError;optinal...)
Computes the hit probability of an air bursting munition.
Returns 3 quantities: 
* The probability of hitting the target defined by'target' with 'p %' number of fragments. 
(p is an optional input, default value is 10). 
* A vector containing the detonation points coordinates
* A vector containing the number of impacting fragments per detonation

The computation are done using MonteCarlo simulations. 
The number of Monte Carlo simulations is N. It is an optional input its default value is 100.
Inputs are:
* tank is a tank object defined using WeaponSystems package
* t is a target object (ExternalBallistics) it contains the position (range of the target)
* proj is a projectile object which contains the projectile characteristics (ExternalBallistics)
* w is the wind magnitude (cross and range wind)
* atmosphere contains the air characteristics
* fuze defines the number of revolutions before detonating the projectile
* aero contains the aerodynamic coefficient of the projectile
* zones is a list of zones that defines the properties of the generated fragments  
* shapes defines the shapes of the fragments 
* dfError is a dictionary which contains the errors on input parameters that need to be propagated
* linError is a dictionary which contains the errors measured that are tabulated
"""
function phitAbm(tank::Tank,t::AbstractTarget, proj::AbstractPenetrator, w::Wind,atmosphere::Air, fuze::Int64, 
    aero::DataFrame,zones::ZoneList, shapes::Array{FragShapes,1},dfError::Dict{String, Float64},linError::Dict{String, Error};N=100,p=10)

weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
weapon.lw = tank.canon.lw
weapon.X2w = tank.rWY
weapon.QE = tank.canon.θ #updates the canon
weapon.AZ = tank.turret.ξ
d1 = Normal(tank.hull.Φ,dfError["cant"])
    
x1 = rand(d1, N)

d2 = Normal(0.0,dfError["cant_measurement"])

x2 = rand(d2, N)

dcross = Normal(w.cross,dfError["crossWind"])
xcross = rand(dcross, N)

drange = Normal(w.range,dfError["rangeWind"])
xrange = rand(drange, N)

dtarget = Normal(t.ρ,dfError["targetRange"])
xtarget = rand(dtarget, N)

dvel = Normal(weapon.u₀,dfError["muzzleVelocity"])
xvel = rand(dvel, N)

dtemp = Normal(atmosphere.t,dfError["airTemperature"])
xtemp = rand(dtemp, N)

abm = createProjectile(proj.mass,proj.calibre)

hits = []
nfrs = 0.0
det = []
RLos = angle_to_dcm(0, -deg2rad(tank.sight.θ), deg2rad(tank.sight.ξ), :XZY)
    
for i=1:N
     cant = x1[i]+x2[i] 
     crossW = xcross[i]
     rangeW=xrange[i]
     ρtarget = xtarget[i]
     mvel = xvel[i] 
     temp = xtemp[i]
     tank.hull.Φ = cant
     w.cross = crossW
     w.range = rangeW
     t.ρ = ρtarget
     t.position = targetPos(t, tank)
     weapon.u₀ = mvel
     proj.velocity=muzzleVel(tank)
     proj.position=muzzlePos(tank)
     #atmosphere.t = temp
     windIn = wind(tank, w)
    # impactP, impactV, tof, αₑ,spin,rounds= trajectoryMPMM(proj, t, weapon,aero)
    variable_bias = [0.0,rand(Normal(0.0,linError["jump"].vertical)),rand(Normal(0.0,linError["jump"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["fireControl"].vertical)),rand(Normal(0.0,linError["fireControl"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["boresight"].vertical)),rand(Normal(0.0,linError["boresight"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["optical"].vertical)),rand(Normal(0.0,linError["optical"].horizontal))]
    
    radom_err =[0.0,rand(Normal(0.0,linError["dispers"].vertical)),rand(Normal(0.0,linError["dispers"].horizontal))] 
 
    detP, detV, detTof, detαₑ,detSpin,detRounds= trajectoryMPMM(proj, t, weapon,aero,tf =fuze,w_bar =windIn)
    #println(i, "\t", ρtarget, "\t", norm(impactV), "\t", tof, "\t", rounds/(2*pi), "\t", mvel)
    detP = inv(RLos)*detP
 
    abm.position = detP+[0.0,linError["fixedBias"].vertical,linError["fixedBias"].horizontal] +variable_bias + radom_err
    abm.position = RLos* abm.position

    abm.velocity = detV
    abm.αₑ = detαₑ
    abm.spin = detSpin

    fragments = frag(abm, zones,shapes)
 
    frs,rays = raytracing(fragments)
 
    #tVulnerable = targetPosition(tVulnerable,t.position)#Vulnerability.target(t,0.5,0.5,0.5)
    tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
 
    hit,shots = shotlines(tVulnerable,rays)
    nfrs = length(frs)
 
    #println(i, "\t", hit, "\t", length(frs))
    push!(hits,hit)
    push!(det,abm.position)
   
    #toto = Dict(zip(rays,frs))
 
    #get(toto,shots[1],"no hit")
 
 
 
 
    #@show tof
 end
 
 phit = length(hits[hits .>= (nfrs*p/100)])/N

 return phit,det,hits
 
end 

"""
sshpAbm(tank,t, proj, fuze, aero,zones, shapes,ϵ;optional....)
Computes the hit probability of an air bursting munition.
Returns 3 quantities: 
* The probability of hitting the target defined by'target' with 'p %' number of fragments. 
(p is an optional input, default value is 10). 
* A vector containing the detonation points coordinates
* A vector containing the number of impacting fragments per detonation

The computation are done using MonteCarlo simulations. 
The number of Monte Carlo simulations is N. It is an optional input its default value is 100.
Inputs are:
* tank is a tank object defined using WeaponSystems package
* t is a target object (ExternalBallistics) it contains the position (range of the target)
* proj is a projectile object which contains the projectile characteristics (ExternalBallistics)
* fuze defines the number of revolutions before detonating the projectile
* aero contains the aerodynamic coefficient of the projectile
* zones is a list of zones that defines the properties of the generated fragments  
* shapes defines the shapes of the fragments 
* ϵ contains the total error defined using a 'SpheError' object.
"""
function sshpAbm(tank::Tank,t::AbstractTarget, proj::AbstractPenetrator, fuze::Float64, targetSize::Vector{Float64},
    aero::DataFrame,zones::ZoneList, shapes::Array{FragShapes,1},ϵ::SpheError;N=100,p=10)
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    
#for i=1:N
RLos = angle_to_dcm(0, -deg2rad(tank.sight.θ), deg2rad(tank.sight.ξ), :XZY)
abm = createProjectile(proj.mass,proj.calibre)
hits = []
nfrs = 0.0
det = []
#tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
tVulnerable = Vulnerability.target(t,targetSize[1],targetSize[2],targetSize[3])
 
    detP, detV, detTof, detαₑ,detSpin,detRounds= trajectoryMPMM(proj, t, weapon,aero,tf =fuze)
    #println(i, "\t", ρtarget, "\t", norm(impactV), "\t", tof, "\t", rounds/(2*pi), "\t", mvel)
    detP = inv(RLos)*detP
 for i=1:N
    abm.position = detP+[ϵ.μ_x,ϵ.μ_y,ϵ.μ_z] +[rand(Normal(0.0,ϵ.σ_x)),rand(Normal(0.0,ϵ.σ_y)), rand(Normal(0.0,ϵ.σ_z))]
    abm.position = RLos* abm.position

    abm.velocity = detV
    abm.αₑ = detαₑ
    abm.spin = detSpin
    fragments = frag(abm, zones,shapes)
 
    frs,rays = raytracing(fragments)
 
    #tVulnerable = targetPosition(tVulnerable,t.position)#Vulnerability.target(t,0.5,0.5,0.5)
    
 
    hit,shots = shotlines(tVulnerable,rays)
    nfrs = length(frs)
 
    #println(i, "\t", hit, "\t", length(frs))
    push!(hits,hit)
    push!(det,abm.position)
   
    
 end
 
 phit = length(hits[hits .>= (nfrs*p/100)])/N

 return phit,det,hits
 
end 

"""
phit(tank,t, proj, w,atmosphere, fuze, aero,zones, shapes,dfError,linError;optinal...)
Computes the distribution of impact points of an inert projectile.
Returns 10 quantities: 
* The gaussian distributions of the impact points coordinates (3)
* The Gaussian distributions of the impact velocity components (3)
* The Gaussian distributions of the impact yaw of repose components (3)
* The Gaussian distribution of the impact spin rate

The computation are done using MonteCarlo simulations. 
The number of Monte Carlo simulations is N. It is an optional input its default value is 100.
Inputs are:
* tank is a tank object defined using WeaponSystems package
* t is a target object (ExternalBallistics) it contains the position (range of the target)
* proj is a projectile object which contains the projectile characteristics (ExternalBallistics)
* w is the wind magnitude (cross and range wind)
* atmosphere contains the air characteristics
* fuze defines the number of revolutions before detonating the projectile
* aero contains the aerodynamic coefficient of the projectile
* zones is a list of zones that defines the properties of the generated fragments  
* shapes defines the shapes of the fragments 
* dfError is a dictionary which contains the errors on input parameters that need to be propagated
* linError is a dictionary which contains the errors measured that are tabulated
"""

function phit(tank::Tank,t::AbstractTarget, proj::AbstractPenetrator,w::Wind,atmosphere::Air, fuze::Int64, 
    aero::DataFrame,dfError::Dict{String, Float64},linError::Dict{String, Error};N=100)
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
d1 = Normal(tank.hull.Φ,dfError["cant"])
    
x1 = rand(d1, N)

d2 = Normal(0.0,dfError["cant_measurement"])

x2 = rand(d2, N)

dcross = Normal(w.cross,dfError["crossWind"])
xcross = rand(dcross, N)

drange = Normal(w.range,dfError["rangeWind"])
xrange = rand(drange, N)

dtarget = Normal(t.ρ,dfError["targetRange"])
xtarget = rand(dtarget, N)

dvel = Normal(weapon.u₀,dfError["muzzleVelocity"])
xvel = rand(dvel, N)

dtemp = Normal(atmosphere.t,dfError["airTemperature"])
xtemp = rand(dtemp, N)

abm = createProjectile(proj.mass,proj.calibre)


nfrs = 0.0
det = []
impactx = zeros(N)
impacty = zeros(N)
impactz = zeros(N)

velx = zeros(N)
vely = zeros(N)
velz = zeros(N)

yawOfreposex = zeros(N)
yawOfreposey = zeros(N)
yawOfreposez = zeros(N)

spin = zeros(N)
RLos = angle_to_dcm(0, -deg2rad(tank.sight.θ), deg2rad(tank.sight.ξ), :XZY)
    
for i=1:N
     cant = x1[i]+x2[i] 
     crossW = xcross[i]
     rangeW=xrange[i]
     ρtarget = xtarget[i]
     mvel = xvel[i] 
     temp = xtemp[i]
     tank.hull.Φ = cant
     w.cross = crossW
     w.range = rangeW
     t.ρ = ρtarget
     t.position = targetPos(t, tank)
     weapon.u₀ = mvel
     proj.velocity=muzzleVel(tank)
     proj.position=muzzlePos(tank)
     #atmosphere.t = temp
     windIn = wind(tank, w)
    # impactP, impactV, tof, αₑ,spin,rounds= trajectoryMPMM(proj, t, weapon,aero)
    variable_bias = [0.0,rand(Normal(0.0,linError["jump"].vertical)),rand(Normal(0.0,linError["jump"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["fireControl"].vertical)),rand(Normal(0.0,linError["fireControl"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["boresight"].vertical)),rand(Normal(0.0,linError["boresight"].horizontal))]+
    [0.0,rand(Normal(0.0,linError["optical"].vertical)),rand(Normal(0.0,linError["optical"].horizontal))]
    
    radom_err =[0.0,rand(Normal(0.0,linError["dispers"].vertical)),rand(Normal(0.0,linError["dispers"].horizontal))] 
 
    detP, detV, detTof, detαₑ,detSpin,detRounds= trajectoryMPMM(proj, t, weapon,aero,tf =fuze,w_bar =windIn)
    #println(i, "\t", ρtarget, "\t", norm(impactV), "\t", tof, "\t", rounds/(2*pi), "\t", mvel)
    detP = inv(RLos)*detP
    detV = inv(RLos)*detV

    abm.position = detP+[0.0,linError["fixedBias"].vertical,linError["fixedBias"].horizontal] +variable_bias + radom_err
    abm.velocity = detV
    abm.αₑ = detαₑ
    abm.spin = detSpin
    
    impactx[i] = abm.position[1]
    impacty[i] = abm.position[2]
    impactz[i] = abm.position[3]
    
    velx[i] = abm.velocity[1] 
    vely[i] = abm.velocity[2] 
    velz[i] = abm.velocity[3] 
    
    yawOfreposex[i] = abm.αₑ[1]
    yawOfreposey[i] = abm.αₑ[2]
    yawOfreposez[i] = abm.αₑ[3]
    
    spin[i] = abm.spin

 end
 
 distImpactx =fit(Normal, impactx)
 distImpacty = fit(Normal, impacty)
 distImpactz = fit(Normal, impactz)

 distVelx = fit(Normal,velx)
 distVely = fit(Normal,vely)
 distVelz = fit(Normal,velz)

 distYawOfReposex = fit(Normal, yawOfreposex)
 distYawOfReposey = fit(Normal, yawOfreposey)
 distYawOfReposez = fit(Normal, yawOfreposez)

 distSpin = fit(Normal, spin)

 return distImpactx, distImpacty, distImpactz, distVelx, distVely, distVelz, distYawOfReposex, distYawOfReposey, distYawOfReposez, distSpin
 
end 

"""
dispersion(target,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError;optional....)
Computes the fixed bias, variable bias and random errors.
Returns 3 quantities: 
* fixed bias an object of type 'Error' which contains the fixed bias in a vertical plane to the LOS
* Variable bias an object of type 'Error3D' which contains the variable bias. 
The error is expressed in the range direction and two other directions perpendicular to the LOS
* random error object of type 'Error' which contains the random error in a vertical plane perpendicular to the LOS

Some error are obtained from interpolations on measured data others by Monte Carlo propagation.
The number of Monte Carlo simulations is N. It is an optional input its default value is 100.
Each error source is propagated individualy the total error is obtained by assumption of non-correlated errors.
Inputs are:
* target is a target object (ExternalBallistics) it contains the position (range of the target)
* tank is a tank object defined using WeaponSystems package
* proj is a projectile object which contains the projectile characteristics (ExternalBallistics)
* aero contains the aerodynamic coefficient of the projectile
* fixedError is a dictionary which contains the fixed bias in mils
* varibleErrorIn is a dictionary which contains the standard deviation on measured values
* variableErrorOut is a dictionary which contains variable bias in mils measured on plane vertical to the LOS 
* randomError is a dictionary which contains the random error in mils 
Other optional arguments are:
* The wind magnitude, its default value is null.
* The atmosphere default are pressure 101325Pa, temperature =15°C and relative humidity = 80.

"""
function dispersion(target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,variableErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0),N=100,atmosphere = Air(101325,15,80))
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    #w = Wind(0.0,0.0)
    σfbh = fixedError["horizontal"](target.ρ) *1e-3*0.981719*target.ρ#interTable(fixed_bias_range,fixed_bias_horizontal,target.ρ)
    σfbv = fixedError["vertical"](target.ρ)*1e-3*0.981719*target.ρ#interTable(fixed_bias_range,fixed_bias_vertical,target.ρ)

    σcantR,σcantV, σcantH=cantError(variableErrorIn["cant"], target,proj,weapon,tank,aero,N=N)[1:3]
    println("cant error", " ", σcantV, " ", σcantH)

    σcrossₜ,σcrossV,σcrossR= windCrossError(variableErrorIn["crossWind"],w,tank,weapon, proj, target, aero,N=N)[1:3]#.σ #m
    println("cross wind error", " ", σcrossₜ.σ)
    σrangeₜ,σrangeH, σrangeR=windRangeError(variableErrorIn["rangeWind"],w,tank,weapon, proj, target, aero,N=N)[1:3]#.σ #m
    println("range wind error", " ", σrangeₜ.σ)
    

    σjh=variableErrorOut["jumpH"](target.ρ)*1e-3*0.981719*target.ρ#interTable(rangeJ,jumpH,target.ρ)*1e-3*0.981719*target.ρ
    σjv=variableErrorOut["jumpV"](target.ρ)*1e-3*0.981719*target.ρ#interTable(rangeJ,jumpV,target.ρ)*1e-3*0.981719*target.ρ

    σfch=variableErrorOut["FCH"](target.ρ)*1e-3*0.981719*target.ρ#interTable(rangeFC,fcH,target.ρ)*1e-3*0.981719*target.ρ
    σfcv=variableErrorOut["FCV"](target.ρ)*1e-3*0.981719*target.ρ#interTable(rangeFC,fcV,target.ρ)*1e-3*0.981719*target.ρ

    σboresightH= variableErrorOut["boresightH"](target.ρ)*1e-3*0.981719*target.ρ#boresight_retention_horizontal*1e-3*0.981719*target.ρ
    σboresightV= variableErrorOut["boresightV"](target.ρ)*1e-3*0.981719*target.ρ#boresight_retention_vertical*1e-3*0.981719*target.ρ

    #dQE,dAZ=rangeError(target,range_std,tank, weapon, proj, aero)
    σrR,σrV, σrH=rangeError(target,variableErrorIn["range"],tank, weapon, proj, aero,N=N) #parallax compensation not done
    #println("range error computed")

    #σrEl,σvEl, σhEl = elevationError(dQE.σ, target,proj,weapon,tank,aero)[1:3]
    #σrAz,σvAz, σhAz = azimuthError(dAZ.σ, target,proj,weapon,tank,aero)[1:3]
    println("range error", " ", σrR, " ", σrV, " ", σrH) 

    σMVr, σMVvₜ, σMVhₜ = muzzleVelError(variableErrorIn["muzzleVelocity"], target,proj,weapon,tank,aero,N=N)[1:3]
    println("muzzle velocity error", " ", σMVvₜ, " ", σMVhₜ)

    σtempR,σtempV, σtempH = temperatureError(variableErrorIn["temperature"], target,proj,weapon,tank,aero,atmosphere,N=N)[1:3]
    println("temperature error", " ", σtempV, " ", σtempH)

    #σopth = optical_path_bending
    #σoptv = optical_path_bending

    #no information on lay error

    σdisph = randomError["horizontal"](target.ρ)*1e-3*0.981719*target.ρ#interTable(range_dispersion,dispersion_horizontal,target.ρ)
    σdispv = randomError["vertical"](target.ρ)*1e-3*0.981719*target.ρ#interTable(range_dispersion,dispersion_vertical,target.ρ)

    fixed_bias = Error(σfbh,σfbv)
    
    σjump=Error(σjh,σjv)
    σfc = Error(σfch,σfcv)
    σboresight=Error(σboresightH,σboresightV)
    #σopt=Error(σopth,σoptv)

    #variable_bias = variableBias(target,tank, weapon, proj, aero, w,atmosphere,σrange=range_std, σcross_wind=cross_wind_std,
    #σrange_wind=range_wind_std, σjump=σjump, σfc=σfc,σboresight=σboresight,
    #σMv=muzzle_velocity_variation, σtemp=air_temperature_std, σopt=σopt, σcant =cant_std, σcant_meas=cant_meas_std)
    #variable_bias.horizontal = variable_bias.horizontal* 1e-3*0.981719*target.ρ
    #variable_bias.vertical = variable_bias.vertical * 1e-3*0.981719*target.ρ

    variable_bias = Error3D(sqrt(σcrossR.σ^2+σrangeR.σ^2+σrR.σ^2+σMVr.σ^2+σtempR.σ^2+σcantR.σ^2),sqrt(σcantH.σ^2+σcrossₜ.σ^2+σrangeH.σ^2+σrH.σ^2+σMVhₜ.σ^2+σjh^2+σfch^2+σboresightH^2), 
    sqrt(σcantV.σ^2+σcrossV.σ^2+ σrangeₜ.σ^2+σrV.σ^2+σMVvₜ.σ^2+σjv^2+σfcv^2+σboresightV^2))

    random_err = Error(σdisph,σdispv)
    return fixed_bias, variable_bias, random_err

end 