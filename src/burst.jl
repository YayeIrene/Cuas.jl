
function phitBurst!(hit::Vector{Float64})

    for j=2:length(hit)
        #global hit
    
        cant = rand(d1)+ rand(d2)
        crossW = rand(dcross)
        rangeW=rand(drange)
        ρtarget = rand(dtarget)
        mvel = rand(dvel) 
        temp = rand(dtemp)

        tank.hull.Φ = cant
        w.cross = crossW
        w.range = rangeW

        t.position = targetPos(t, tank)
    
        proj.velocity=muzzleVel(tank)
        proj.position=muzzlePos(tank)

        windIn = wind(tank, w)

    variable_bias = [0.0,rand(Normal(0.0,σjump.vertical)),rand(Normal(0.0,σjump.horizontal))]+[0.0,rand(Normal(0.0,σfc.vertical)),rand(Normal(0.0,σfc.horizontal))]+
    [0.0,rand(Normal(0.0,σboresight.vertical)),rand(Normal(0.0,σboresight.horizontal))]+[0.0,rand(Normal(0.0,σopt.vertical)),rand(Normal(0.0,σopt.horizontal))]
    
    radom_err =[0.0,rand(Normal(0.0,σdisph*mils2m)),rand(Normal(0.0,σdisph*mils2m))] 

        detP, detV, detTof, detαₑ,detSpin,detRounds= trajectoryMPMM(proj, t, weapon,aero,tf =2000,w_bar =windIn)
        abm.position = detP+[0.0,fixed_bias.vertical,fixed_bias.horizontal] +variable_bias + radom_err
        abm.velocity = detV
        abm.αₑ = detαₑ
        abm.spin = detSpin
        fragments = frag(abm, [z1,z2],[s1,s2,s3,s4,s5])

        frs,rays = raytracing(fragments)

        tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)

        hit[j],shots = shotlines(tVulnerable,rays)
    end 
    #return hit
end 


function phitAbmBurst(tank::Tank,t::AbstractTarget, weapon::Gun, proj::AbstractPenetrator, tVulnerable::CSGGenerator{Float64}, 
    w::Wind,atmosphere::Air, fuze::Int64, burstLength::Int64,
    aero::DataFrame,zones::Array{Zone,1}, shapes::Array{FragShapes,1},dfError::Dict{String, Float64},linError::Dict{String, Error};N=100,p=10)

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
monteCarloMatrix= Matrix{Float64}(undef,N,burstLength)

nfrs = 0.0
hitsPerBurst = Vector{Float64}(undef,burstLength) 
pnhit = Vector{Float64}(undef,N) 

    
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
     #t.ρ = ρtarget
     t.position = targetPos(t, tank)
     #weapon.u₀ = mvel
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
 
    abm.position = detP+[0.0,linError["fixedBias"].vertical,linError["fixedBias"].horizontal] +variable_bias + radom_err
    abm.velocity = detV
    abm.αₑ = detαₑ
    abm.spin = detSpin
    fragments = frag(abm, zones,shapes)
 
    frs,rays = raytracing(fragments)
 
    #tVulnerable = targetPosition(tVulnerable,t.position)#Vulnerability.target(t,0.5,0.5,0.5)
    tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
 
    hitsPerBurst[1],shots = shotlines(tVulnerable,rays)
    nfrs = length(frs)

    #phitBurst!(hitsPerBurst)

    for j=2:length(burstLength)
        #global hit
    
        cant = rand(d1)+ rand(d2)
        crossW = rand(dcross)
        rangeW=rand(drange)
        ρtarget = rand(dtarget)
        mvel = rand(dvel) 
        temp = rand(dtemp)

        tank.hull.Φ = cant
        w.cross = crossW
        w.range = rangeW

        t.position = targetPos(t, tank)
    
        proj.velocity=muzzleVel(tank)
        proj.position=muzzlePos(tank)

        windIn = wind(tank, w)

    variable_bias = [0.0,rand(Normal(0.0,σjump.vertical)),rand(Normal(0.0,σjump.horizontal))]+[0.0,rand(Normal(0.0,σfc.vertical)),rand(Normal(0.0,σfc.horizontal))]+
    [0.0,rand(Normal(0.0,σboresight.vertical)),rand(Normal(0.0,σboresight.horizontal))]+[0.0,rand(Normal(0.0,σopt.vertical)),rand(Normal(0.0,σopt.horizontal))]
    
    radom_err =[0.0,rand(Normal(0.0,σdisph*mils2m)),rand(Normal(0.0,σdisph*mils2m))] 

        detP, detV, detTof, detαₑ,detSpin,detRounds= trajectoryMPMM(proj, t, weapon,aero,tf =2000,w_bar =windIn)
        abm.position = detP+[0.0,fixed_bias.vertical,fixed_bias.horizontal] +variable_bias + radom_err
        abm.velocity = detV
        abm.αₑ = detαₑ
        abm.spin = detSpin
        fragments = frag(abm, [z1,z2],[s1,s2,s3,s4,s5])

        frs,rays = raytracing(fragments)

        tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)

        hitsPerBurst[j],shots = shotlines(tVulnerable,rays)
    end 
 
    monteCarloMatrix[i,:] = hitsPerBurst

    pnhit[i] = length(hitsPerBurst[hitsPerBurst .>= (nfrs*p/100)])/burstLength

    
 end
 
 phit = mean(pnhit)
 
 #return phit
 
end 

function phitBurst(n::Int64,k::Int64,sshp::Float64) #probability of getting exactly k hits in a burst of n rounds.
    #p=1-(binomial(n,k)*sshp^n*(1-sshp)^(n-k))
    d = Binomial(n,sshp)
    p = pdf(d,k)
end

function phitSingleBurst(n::Int64,sshp::Float64)
    p = 1-phitBurst(n,0,sshp)
end 

function phitAtLeastBurst(n::Int64,l::Int64,sshp::Float64)
    pi =0.0
    for k=0:(l-1)
        pi=pi+phitburst(n,k,sshp)
    end 
    p=1-pi
end 
"""
phitAtLeastOneBurst(n,sshp)

returns the burst hit probability
* n is the number of rounds in the burst
* sshp is the single round hit probability  
"""
function phitAtLeastOneBurst(n::Int64,sshp::Float64)
    p = 1-(1-sshp)^n
end 

function numberOfHitBurst(n::Int64,sshp::Float64)
    μ = sshp*n
    temp =0.0
    for k=0:n 
        temp=temp+k^2*phitburst(n,k,sshp)
    end 
    σ = sqrt(μ^2+temp)
    d =Normal(μ,σ)
end 

function phitBurstAbm(tank::Tank,t::AbstractTarget, proj::AbstractPenetrator, fuze::Float64, targetSize::Vector{Float64},
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

    monteCarloMatrix= Matrix{Float64}(undef,N,burstLength)

nfrs = 0.0
hitsPerBurst = Vector{Float64}(undef,burstLength) 
pnhit = Vector{Float64}(undef,N) 

 for i=1:N
    count = 0
    abm.position = detP+[ϵ.μ_x,ϵ.μ_y,ϵ.μ_z] +[rand(Normal(0.0,ϵ.σ_x)),rand(Normal(0.0,ϵ.σ_y)), rand(Normal(0.0,ϵ.σ_z))]
    abm.position = RLos* abm.position

    abm.velocity = detV
    abm.αₑ = detαₑ
    abm.spin = detSpin
    fragments = frag(abm, zones,shapes)
 
    frs,rays = raytracing(fragments)
 
    #tVulnerable = targetPosition(tVulnerable,t.position)#Vulnerability.target(t,0.5,0.5,0.5)
    
 
    #hit,shots = shotlines(tVulnerable,rays)
    hitsPerBurst[1],shots = shotlines(tVulnerable,rays)
    nfrs = length(frs)
 
    #println(i, "\t", hit, "\t", length(frs))
    #push!(hits,hit)
    #push!(det,abm.position)
    for j=2:n 
        abm.position = abm.position
    end
   
    
 end
 
 phit = length(hits[hits .>= (nfrs*p/100)])/N

 return phit,det,hits
 
end 
