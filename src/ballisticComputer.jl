
@resumable function shoot(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, shapes::Array{FragShapes,1}, 
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), deltaR=5, N=100,p=1)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze, " "," rounds")
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero, tf =fuze)

     fixed_bias, variable_bias, random_err=dispersion(t,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
     #sshpAbm(tank,detonation, proj, fuze, aero,zones,ϵ,N=N)[1]
     Pnhit=Pnhit*(1-sshp)
     burst=burst+1
     println("Hit probability is"," ", sshp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
   println("Target missed!")
end
end 
@resumable function shoot(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze, " "," rounds")
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
    
     sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     #burstp = phitAtLeastOneBurst(n,sshp)
     Pnhit=Pnhit*(1-sshp)
     #Pnhit=Pnhit*(1-burstp)
     burst=burst+1
     println("error", " ",ϵ)
     println("sshp", " ", sshp)
     #println("Hit probability is"," ", burstp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    target.alive=false
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
   println("Target missed!")
end
end 


#Multiple
@resumable function shoot(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,targetsList::Vector{Tuple{Float64, Float64}},
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    drones=[]
    Pnhiti = Vector{Float64}(undef,length(targetsList))
    Pnhiti .=1.0
    for tar in targetsList 
        push!(drones,ExternalBallistics.createTarget())
    end
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze)
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
i=1

     for drone in drones
        #sshp = sshpAbm(tank,t,detonation, proj, aero,zd,ϵ,N=100)[1]
        drone.position = t.position+[targetsList[i][1], 0.0, targetsList[i][2]]
       sshp= sshpAbm(tank,drone,detonation, proj, aero,zones,ϵ,N=N)[1]
        Pnhiti[i]=Pnhiti[i]*(1-sshp)
        println("phit drone $i", " ", 1-Pnhiti[i])
        i+=1
        #push!(phits,sshp)
        end 
    
     #sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     #Pnhit=Pnhit*(1-sshp)
     Pnhit = minimum(Pnhiti)
     burst=burst+1
     #println("Hit probability is"," ", sshp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
    missed = count(j->(j<= 0.95), 1-Pnhiti)
   println("$missed Target missed!")
end
end 

#Multiple burst
@resumable function shoot(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,n::Int64,targetsList::Vector{Tuple{Float64, Float64}},
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    drones=[]
    Pnhiti = Vector{Float64}(undef,length(targetsList))
    Pnhiti .=1.0
    for tar in targetsList 
        push!(drones,ExternalBallistics.createTarget())
    end
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze)
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
i=1

     for drone in drones
        #sshp = sshpAbm(tank,t,detonation, proj, aero,zd,ϵ,N=100)[1]
        drone.position = t.position+[targetsList[i][1], 0.0, targetsList[i][2]]
       sshp= sshpAbm(tank,drone,detonation, proj, aero,zones,ϵ,N=N)[1]
       burstp = phitAtLeastOneBurst(n,sshp)
     #Pnhit=Pnhit*(1-sshp)
     Pnhiti[i]=Pnhiti[i]*(1-burstp)
       # Pnhiti[i]=Pnhiti[i]*(1-sshp)
        println("phit drone $i", " ", 1-Pnhiti[i])
        i+=1
        #push!(phits,sshp)
        end 
    
     #sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     #Pnhit=Pnhit*(1-sshp)
     Pnhit = minimum(Pnhiti)
     burst=burst+1
     #println("Hit probability is"," ", sshp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
    missed = count(j->(j<= 0.95), 1-Pnhiti)
   println("$missed Target missed!")
end
end 


#Burst
@resumable function shoot(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,n::Int64,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze, " "," rounds")
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
    
     sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     burstp = phitAtLeastOneBurst(n,sshp)
     #Pnhit=Pnhit*(1-sshp)
     Pnhit=Pnhit*(1-burstp)
     burst=burst+1
     println("error", " ",ϵ)
     println("sshp", " ", sshp)
     println("Hit probability is"," ", burstp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    target.alive=false
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
   println("Target missed!")
end
end 
#-----------------------------------------------------------------------------------------------------------------------

#Multiple burst Mobile/Mobile
@resumable function shoot(sim::Simulation,target::Uav, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,n::Int64,targetsList::Vector{Tuple{Float64, Float64}},
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    drones=[]
    Pnhiti = Vector{Float64}(undef,length(targetsList))
    Pnhiti .=1.0
    for tar in targetsList 
        push!(drones,ExternalBallistics.createTarget())
    end
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     #t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle

     leadAim = lead(target, tof) # computes the lead aimpoint
     leadAimSph  = SphericalFromCartesian()([leadAim[1],-leadAim[3],leadAim[2]])
     tank.sight.θ = rad2deg(leadAimSph.ϕ)
     tank.sight.ξ= rad2deg(leadAimSph.θ)
     t.ρ = leadAimSph.r

     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze)
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
i=1

     for drone in drones
        #sshp = sshpAbm(tank,t,detonation, proj, aero,zd,ϵ,N=100)[1]
        drone.position = t.position+[targetsList[i][1], 0.0, targetsList[i][2]]
       sshp= sshpAbm(tank,drone,detonation, proj, aero,zones,ϵ,N=N)[1]
       burstp = phitAtLeastOneBurst(n,sshp)
     #Pnhit=Pnhit*(1-sshp)
     Pnhiti[i]=Pnhiti[i]*(1-burstp)
       # Pnhiti[i]=Pnhiti[i]*(1-sshp)
        println("phit drone $i", " ", 1-Pnhiti[i])
        i+=1
        #push!(phits,sshp)
        end 
    
     #sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     #Pnhit=Pnhit*(1-sshp)
     Pnhit = minimum(Pnhiti)
     burst=burst+1
     #println("Hit probability is"," ", sshp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    target.alive=false
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
    missed = count(j->(j<= 0.95), 1-Pnhiti)
   println("$missed Target missed!")
end
end 

#Mobile/mobile
@resumable function shoot(sim::Simulation,target::Uav, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    Δ_fireControl::Float64, Δ_fireDemand::Float64, Δ_fireDelay::Float64, Δ_shotExit::Float64,
    zones::ZoneList, fuze::Fuze,n::Int64,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict;w = Wind(0.0,0.0), N=100)
    Pnhit=1.0
    burst=0
    t = ExternalBallistics.createTarget() #creates ExternalBallistics target
    weapon = createGun(tank.canon.u₀,tank.latitude,tank.canon.θ,tank.turret.ξ,tank.canon.twist)
    weapon.lw = tank.canon.lw
    weapon.X2w = tank.rWY
    weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
    # t.position = targetPos(t, tank)
   while target.ρ > 0 && (1-Pnhit) < 0.95
    @yield timeout(sim, Δ_fireControl)
     println("Fire Control solution for target"," ", now(sim))
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     weapon.QE = tank.canon.θ
     weapon.AZ = tank.turret.ξ
     t.ρ = target.ρ
     t.position = targetPos(t, tank) 
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank)
     #println("ballistic correction", " ", QE, " ", AZ)
      #println("tank ballistic correction", " ", tank.canon.θ, " ", tank.turret.ξ)
      #println("weapon ballistic correction", " ", weapon.QE, " ", weapon.AZ)

     #tank.canon.θ = QE
     #tank.turret.ξ = AZ

     #weapon.QE = QE
     #weapon.AZ = AZ
     proj.position=muzzlePos(tank)
     proj.velocity=muzzleVel(tank)
     tof,αₑ,spin,rounds = trajectoryMPMM(proj, t, weapon,aero)[3:6]

     #t.ρ,  tank.sight.θ, tank.sight.ξ  = lead(target, tof) # computes the lead angle
     leadAim = lead(target, tof) # computes the lead aimpoint
     leadAimSph  = SphericalFromCartesian()([leadAim[1],-leadAim[3],leadAim[2]])
     tank.sight.θ = rad2deg(leadAimSph.ϕ)
     tank.sight.ξ= rad2deg(leadAimSph.θ)
     t.ρ = leadAimSph.r

     tank.canon.θ = tank.sight.θ #moves the canon to the lead angle
     tank.turret.ξ = tank.sight.ξ 
     #target.position = targetPos(target, tank) # computes the new target position
     #t.ρ = target.ρ
     t.position = targetPos(t, tank)

     proj.position=muzzlePos(tank) #updates muzzle position
     proj.velocity=muzzleVel(tank) #updates muzzle velocity
     weapon.QE = tank.canon.θ #updates the canon
     weapon.AZ = tank.turret.ξ
     #QE,AZ=QEfinderMPMM!(target, proj, weapon,aero)
     adjustedFire!(t, proj, weapon,aero,tank) # adjust for the lead angle

     println("lead angle", " ", weapon.QE, " ", weapon.AZ)

     #QE,AZ=QEfinderMPMM!(t, proj, weapon,aero)    



     @yield timeout(sim, Δ_fireDemand)
     println("Fire at target !"," ", now(sim))

     @yield timeout(sim, Δ_fireDelay)
     println("firing delay for target"," ", now(sim))

     @yield timeout(sim, Δ_shotExit)
     println("Shot exit"," ", now(sim))
     if target.ρ <10
        break 
     end 
     
     #tf = rounds-20# 0.048995 #timefuze


     #impactP, impactV, tof, αₑ,spin = trajectoryMPMM(proj, target, weapon,aero, tspan=(0,tof-tf))
     #fuze =  rounds/(2*pi)-deltaR#5#10#2000
     println("fuze set to ", " ", fuze, " "," rounds")
     detonation = deepcopy(t)
    detonation.position = t.position - [fuze.offset,fuze.height,0.0]
     impactP, impactV, tof, αₑ,spin,rounds = trajectoryMPMM(proj, detonation, weapon,aero)

     fixed_bias, variable_bias, random_err=dispersion(detonation,tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,N=N)

      #target.ρ= track.ρ
     ϵ = SpheError(variable_bias.range,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #ϵ = SpheError(0.0,sqrt(variable_bias.vertical^2+random_err.vertical^2),sqrt(variable_bias.horizontal^2+random_err.horizontal^2),0.0,fixed_bias.vertical,fixed_bias.horizontal)
     #sshp = SSHP(target,ϵ)
     #tVulnerable = Vulnerability.target(t,0.5,0.5,0.5)
     #sshp = sshpAbm(tank,t, proj, fuze,target.size,aero,zones,shapes,ϵ,p=p,N=N)[1]
    
     sshp = sshpAbm(tank,t,detonation, proj, aero,zones,ϵ,N=N)[1]
     burstp = phitAtLeastOneBurst(n,sshp)
     #Pnhit=Pnhit*(1-sshp)
     Pnhit=Pnhit*(1-burstp)
     burst=burst+1
     println("error", " ",ϵ)
     println("sshp", " ", sshp)
     println("Hit probability is"," ", burstp)

     @yield timeout(sim, tof)
     println("Target strike"," ", now(sim))
     #targetPosition!(target,now(sim))
    
end
Phit = (1-Pnhit)
if target.ρ >0
    target.alive=false
    println("Target position at hit"," ", target.ρ)
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
else
    println("Hit probability is"," ", Phit)
    println("burst length is", " ", burst)
   println("Target missed!")
end
end 