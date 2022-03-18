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

function lead(target::Uav,Δt::Float64)
    x = target.position[1] + target.velocity[1]*Δt
    y = target.position[2] + target.velocity[2]*Δt
    z = target.position[3] + target.velocity[3]*Δt
    return x,y,z
end


@resumable function fly!(sim::Simulation,target::AbstractTarget,Δt::Float64)

    while target.ρ > 0 && target.alive
         @yield timeout(sim, Δt)
        target.ρ = target.ρ - target.ρ′*Δt
        target.θ = target.θ + target.θ′*Δt
        target.ϕ = target.ϕ + target.ϕ′*Δt
        println("target position ", " ",target.ρ)

    end
end


@resumable function fly!(sim::Simulation,target::Uav,Δt::Float64)

    while target.ρ > 0 && target.alive
         @yield timeout(sim, Δt)
        target.position[1] = target.position[1] + target.velocity[1]*Δt
        target.position[2] = target.position[2] + target.velocity[2]*Δt
        target.position[3] = target.position[3] + target.velocity[3]*Δt
        println("target position ", " ",target.position)

    end
end

@resumable function fly!(sim::Simulation,target::Uav, tank::Tank,Δt::Float64)
    #target.ρ = norm(target.position-tank.position)
    aiming!(target,tank)
    #engagementCart= target.position-tank.position
    #engagementSph  = SphericalFromCartesian()([engagementCart[1],-engagementCart[3],engagementCart[2]])
    #target.θ = rad2deg(90.0-engagementSph.ϕ)
    #target.ϕ = rad2deg(engagementSph.θ)

    #target.ρ = engagementSph.r

    while target.ρ > 0 && target.alive
         @yield timeout(sim, Δt)
        target.position[1] = target.position[1] + target.velocity[1]*Δt
        target.position[2] = target.position[2] + target.velocity[2]*Δt
        target.position[3] = target.position[3] + target.velocity[3]*Δt
        #target.ρ = norm(target.position-tank.position)
        aiming!(target,tank)
        #engagementCart= target.position-tank.position
        #engagementSph  = SphericalFromCartesian()([engagementCart[1],-engagementCart[3],engagementCart[2]])
        #target.θ = rad2deg(90.0-engagementSph.ϕ)
        #target.ϕ = rad2deg(engagementSph.θ)

        #target.ρ = engagementSph.r
        println("target position ", " ",target.ρ)

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
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, shapes,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w,deltaR=deltaR, N=N,p=p)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

@resumable function acquisition(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze;
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
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

#Multiple
@resumable function acquisition(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze, targetsList::Vector{Tuple{Float64, Float64}};
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
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,targetsList,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

#For a burst
@resumable function acquisition(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze, n::Int64;
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
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,n,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

@resumable function advance!(sim::Simulation,tank::Tank, target::AbstractTarget, Δt::Float64)
  
    while target.alive
        @yield timeout(sim, Δt)
        tank.position[1] = tank.position[1] +  tank.velocity[1]*Δt
        tank.position[2] = tank.position[2] +  tank.velocity[2]*Δt
        tank.position[3] = tank.position[3] +  tank.velocity[3]*Δt
       println("tank position ", " ",tank.position)

   end
    
    
end 

function aiming!(target::Uav,tank::Tank)
    los= target.position-tank.position
    aimpoint  = SphericalFromCartesian()([los[1],-los[3],los[2]])
    target.θ = rad2deg(aimpoint.ϕ)
    target.ϕ = rad2deg(aimpoint.θ)
    target.ρ = aimpoint.r
end 

#For a burst Multiple
@resumable function acquisition(sim::Simulation,target::AbstractTarget, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze, n::Int64,  targetsList::Vector{Tuple{Float64, Float64}};
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
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,n,targetsList,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

#--------------------------------------------------------------------------------------------------------------------

#For a burst Multiple Mobile/Mobile
@resumable function acquisition(sim::Simulation,target::Uav, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze, n::Int64,  targetsList::Vector{Tuple{Float64, Float64}};
    Δt_fireCommand = 1.0, Δ_slew = 2.0, Δ_aim = 2.0, Δ_range= 0.0,
    Δ_fireControl = 1.0, Δ_fireDemand = 0.0, Δ_fireDelay = 0.0, Δ_shotExit = 1.0,
    w = Wind(0.0,0.0),deltaR=5, N=100,p=1)

    println("Fire at target ")
  
    println("Target position at fire command"," ", target.ρ)
    println("Tank position at fire command", " ", tank.position)


    @process fly!(sim,target,1.0)
    @process advance!(sim,tank, target, 1.0)

    @yield timeout(sim, Δ_slew)
    println("Slew weapon at target"," ", now(sim))
    #engagementCart= target.position-tank.position
    #engagementSph  = SphericalFromCartesian()([engagementCart[1],-engagementCart[3],engagementCart[2]])
    #target.θ = rad2deg(90.0-engagementSph.ϕ)
    #target.ϕ = rad2deg(engagementSph.θ)
    #target.ρ = engagementSph.ρ


    tank.turret.ξ = target.ϕ
    

    #@yield timeout(sim, Δ_aim)
    #println("Aim at target $id"," ", now(sim))
    tank.canon.θ = target.θ



    @yield timeout(sim, Δ_range)
    println("Range at target"," ", now(sim))
   
    tank.sight.θ = target.θ
    tank.sight.ξ = target.ϕ
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,n,targetsList,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end

#Mobile/mobile
@resumable function acquisition(sim::Simulation,target::Uav, tank::Tank, proj::AbstractPenetrator, aero::DataFrame,
    fixedError::Dict,varibleErrorIn::Dict,variableErrorOut::Dict,randomError::Dict,
    zones::ZoneList, fuze::Fuze, n::Int64;
    Δt_fireCommand = 1.0, Δ_slew = 2.0, Δ_aim = 2.0, Δ_range= 0.0,
    Δ_fireControl = 1.0, Δ_fireDemand = 0.0, Δ_fireDelay = 0.0, Δ_shotExit = 1.0,
    w = Wind(0.0,0.0),deltaR=5, N=100,p=1)

    
    println("Fire at target ")
  
    println("Target position at fire command"," ", target.ρ)
    println("Tank position at fire command", " ", tank.position)


    @process fly!(sim,target,1.0)
    @process advance!(sim,tank, target, 1.0)

    @yield timeout(sim, Δ_slew)
    println("Slew weapon at target"," ", now(sim))
    #engagementCart= target.position-tank.position
    #engagementSph  = SphericalFromCartesian()([engagementCart[1],-engagementCart[3],engagementCart[2]])
    #target.θ = rad2deg(90.0-engagementSph.ϕ)
    #target.ϕ = rad2deg(engagementSph.θ)
    #target.ρ = engagementSph.ρ


    tank.turret.ξ = target.ϕ
    

    #@yield timeout(sim, Δ_aim)
    #println("Aim at target $id"," ", now(sim))
    tank.canon.θ = target.θ



    @yield timeout(sim, Δ_range)
    println("Range at target"," ", now(sim))
   
    tank.sight.θ = target.θ
    tank.sight.ξ = target.ϕ
    

    @yield @process shoot(sim,target,tank, proj, aero,Δ_fireControl, Δ_fireDemand, Δ_fireDelay, Δ_shotExit,
    zones, fuze,n,fixedError,varibleErrorIn,variableErrorOut,randomError,w=w, N=N)     

    @yield timeout(sim, 2.0)
    println("Target damage assessment"," ", now(sim))

    @yield timeout(sim, 1.0)
    println("Reengage Target?"," ", now(sim))


    #@yield release(bcs) # customer exits service
    println("Next target : ", now(sim))

end