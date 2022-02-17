function targetPosition!(target::AbstractTarget,t::Float64)

    if target.ρ >0

    Δt = t - target.t0

        target.ρ = target.ρ - target.ρ′*Δt
        #println("target $id position ", " ",target.ρ)
    end


end

function lead(target::AbstractTarget, Δt::Float64)
    #t=deepcopy(target)
    
    ρ = target.ρ - target.ρ′*Δt
    θ = target.θ + target.θ′*Δt
    ϕ = target.ϕ + target.ϕ′*Δt
    #QE,AZ = QEfinderMPMM!(t, proj, weapon,aero)
    return ρ,θ,ϕ
end 

@resumable function fly!(sim::Simulation,target::AbstractTarget,Δt::Float64)

    while target.ρ > 0
         @yield timeout(sim, Δt)
        target.ρ = target.ρ - target.ρ′*Δt
        target.θ = target.θ + target.θ′*Δt
        target.ϕ = target.ϕ + target.ϕ′*Δt
        println("target position ", " ",target.ρ)


    end
end

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
#=
"""
    acquisition(sim,target, tank, proj, aero,fixedError,varibleErrorIn,variableErrorOut,randomError,zones, shapes; optional arguments)

    Is a simJulia process which simulate the firing chain process for an ABM projectile. From target acquisition to damage assessment.
    Input parameters are :
    * sim is the simulation environment
    * target is a target object which contains target information
    * tank is a tank object (WeaponSystems) 
    * proj contains the projectile informaiton
    * aero contains the aerodynamic coefficient
    * fixedError is a dictionary which contains the fixed bias (in mils) measured in a plane perpendicular to the LOS
    * varibleErrorIn is a dictionary which contains the standard deviation of measured quantities that are propagated via Monte Carlo
    * variableErrorOut is a dictionary which contains variable bias (in mils) measured in a plane perpendicular to the LOS 
    * randomError is a dictionary which contains random errors (in mils) measured in a plane perpendicular to the LOS
    * zones is an object of ZoneList (ExternalBallistics) which contains informaition on generated fragments
    * shapes is a vector which contains the shapes of the fragments

    Optional arguments are:
    * Δt_fireCommand is the time delay for issuing the fire command default value is 1.0
    * Δ_slew is the time delay to lay the canon in a desired position. The default value is 2.0
    * Δ_aim is the time delay to aim at the target. The default value is  2.0
    * Δ_range is the time delay to get the target range. The default value is 0
    * Δ_fireControl is the time delay to get the ballistic correction. The default value is 1.0
    * Δ_fireDemand is the time delay for the fire demand. The default value is 0.0
    * Δ_fireDelay is the time delay for the fire. The default value is  0.0
    * Δ_shotExit is the time delay for the shot to exit. The default value is 1.0
    * w is the wind (cross and range). The default value is  Wind(0.0,0.0)
    * deltaR defines the detonation distance. It expressed in projectile revolution. To the number of revolutions achieved 
    to hit the target, deltaR is soustracted. The default value is 5
    * N is the number of Monte Carlo simulations. The defaults is 100
    * p defines the number in % of fragment that hit target to be considered a hit

"""
=#
@resumable function acquisition(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, shapes::Array{FragShapes,1};
    Δt_fireCommand = 1.0, Δ_slew = 2.0, Δ_aim = 2.0, Δ_range= 0.0,
    Δ_fireControl = 1.0, Δ_fireDemand = 0.0, Δ_fireDelay = 0.0, Δ_shotExit = 1.0,
    w = Wind(0.0,0.0),deltaR=5, N=100,p=1)

    #weapon = createGun(u₀,latitude,canon.θ,tourelle.ξ,tc)
    #weapon.lw = canon.lw
    #weapon.X2w = tank.rWY
    #weapon.QE = tank.canon.θ #updates the canon
     #weapon.AZ = tank.turret.ξ

    #@yield timeout(sim, Δt_fireCommand)
    println("Fire at target ")
  
    println("Target position at fire command"," ", target.ρ)


    @process fly!(sim,target,1.0)

    @yield timeout(sim, Δ_slew)
    println("Slew weapon at target"," ", now(sim))
    tank.turret.ξ = target.ϕ
    

    #@yield timeout(sim, Δ_aim)
    #println("Aim at target $id"," ", now(sim))
    tank.canon.θ = target.θ



    @yield timeout(sim, Δ_range)
    println("Range at target"," ", now(sim))
   
    tank.sight.θ = target.θ
    tank.sight.ξ = target.ϕ
    

    @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, shapes,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,deltaR=deltaR, N=N,p=p)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end