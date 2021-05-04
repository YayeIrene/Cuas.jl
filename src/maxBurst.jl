module phit
using .types


    elseif sim_burst == 1
        #VEND = Float64[]
        #sshp, vend = SSHP(velocity_max, target, D)
        #println("toto1")
        p = Projectile(0.0, velocity, mass, calibre)
        flight_time = @yield @process projectile_flight(sim, p, engagement_distance, delta_t, delta_x)
        println("projectile velocity", " ", p.velocity)
        println("projectile position", " ", p.position)
        Std=σtot(p.position, flight_time)
        #r = abs(rand(Truncated(Normal(0,Std),-1,1)))
        r = rand(Normal(0,Std))
        #p.position = p.position + r
        #println("updated projectile position", " ", p.position)

        #if p.position > engagement_distance
        #    println("target missed")
            #break
        #end
        #push!(VEND, p.velocity)
        #v_rot = v_rot_ini/velocity*p.velocity
        #fragment_velocity = sqrt(p.velocity^2 + v_rot^2 )
        vp =perpendicular_velocity(p, velocity)
        fragment_velocity = sqrt(p.velocity^2 + vp^2 )
        #fragment_velocity = fragment_vel(p, velocity)
        f = Fragment(0, fragment_velocity, fragment_mass, fragment_diameter)
        #detonation_param = @yield @process detonation(sim, p, flight_time, target, f, v_rot)
        d = engagement_distance - p.position
        println("detonation distance", " ", d)
        f_flight = @yield @process fragment_flight(sim, f, d, delta_t)
        #println(" ")
        det_param = @yield @process detonation(sim, p, flight_time, target, f, vp, d)
        n_disp = det_param[1]
        r_lethal = det_param[2]
        println("fragment velocity", " ", f.velocity, " ", "fragment prosition", " ", f.position)
        #d = detonation_param[2]
        println("lethal radius", " ", r_lethal)
        size = target.d + r_lethal*2

        Phiti =  SSHP(size, p.position, flight_time)
        println("sshp", " ", Phiti)


        println("fragment velocity", " ", f.velocity)

        Pkfi = kill_probability(f)
        println("pki for fragment", " ", Pkfi )
        #n_disp = detonation_param[1]
        println("number of fragments", " ", n_disp)
        Pki = 1-(1-Pkfi)^n_disp
        Pkhi = Pki*Phiti
        println(" pki for projectile", " ", Pki)
        println(" pkhi for projectile", " ", Pkhi)
        #push!(VEND, vend)
        if Pkhi > HIT_PROBABILITY
            burst = 1
        else
            burst = Int(ceil(log(1-HIT_PROBABILITY)/log(1-Pkhi)))
        end
        if burst > MAX_BURST_SIZE
            burst = MAX_BURST_SIZE
        end
        Pktot = 1-(1-Pkhi)^burst
        #Phit = 1-(1-sshp)^burst
        #Ekinb = burst * 0.5*mass*(velocity_max)^2
        Ekinb = burst * 0.5*mass*(velocity)^2
        #burst, Phit, Ekinb, mean(VEND)
        burst, Pktot, Ekinb
    elseif sim_burst == 2

    #    while (1-Pkhn) < HIT_PROBABILITY && burst<MAX_BURST_SIZE
            p = Projectile(0.0, v_min, mass, calibre)
            #flight_time = @yield @process projectile_flight(sim, p, engagement_distance, 0.01)
            Ttot = @yield @process projectile_flight(sim, p, engagement_distance, delta_t, delta_x)

            #burst = 3
            println(Ttot)
            #size = target.d
            i= 1

            #for i in 1:burst
            while (1-Pkhn) < HIT_PROBABILITY && burst<MAX_BURST_SIZE
            #while (1-Pkn) < HIT_PROBABILITY && burst<MAX_BURST_SIZE
            #    println(i)
            #end
                Tfi = Ttot - (i - 1)/RoF
                println("tfi", " ", Tfi)
                vi = @yield @process intelligent_burst(sim,Tfi, engagement_distance, v_max, Ttot)
                println("velocity: $vi")
                if vi == 0.0
                    println("velocity te hoog")
                    #burst = burst-1
                    #@goto End
                    break
                end
#                sshp, vend = SSHP_lu(Projectile_lu(0.0, velocity, mass, calibre), Target_lu(D, v_target, L_target, D_target, δ), R_let, σball, σatm, Verror)

                pi = Projectile(0.0, vi, mass, calibre)
                tfi = @yield @process projectile_flight(sim, pi, engagement_distance, delta_t, delta_x)
                println("projectile velocity", " ", pi.velocity)
                println("projectile position", " ", pi.position)

                Std=σtot(pi.position, tfi)
                #r = abs(rand(Truncated(Normal(0,Std),-1,1)))
                r = rand(Normal(0,Std))
                pi.position = pi.position + r
                println("updated projectile position", " ", pi.position)
                #v_roti = v_rot_ini/velocity*pi.velocity
                #fragment_velocity = sqrt(pi.velocity^2 + v_roti^2 )
                if pi.position < engagement_distance
                vp =perpendicular_velocity(pi, vi)
                fragment_velocity = sqrt(pi.velocity^2 + vp^2 )
                #fragment_velocity = fragment_vel(p, vi)
                fi = Fragment(0, fragment_velocity, fragment_mass, fragment_diameter)
                d = engagement_distance - pi.position
                println("detonation distance", " ", d)
                det_param = @yield @process detonation(sim, pi, tfi, target, fi, vp, d)
                n_dispi = det_param[1]
                r_lethal = det_param[2]
                println("fragment velocity", " ", fi.velocity, " ", "fragment prosition", " ", fi.position)
                println("lethal radius", " ", r_lethal)
                size = target.d + r_lethal*2

                Phiti =  SSHP(size, pi.position, tfi)
                #println("sshp: $sshp, vend: $vend")
                println("flight time"," ", tfi)
                println("Phit", " ",Phiti )
#                push!(VEND, vend)
                Ekinb = Ekinb + 0.5*mass*(vi)^2
                Pnhit=Pnhit*(1-Phiti)
                 @yield @process fragment_flight(sim, fi, d, 0.01)
                 println("fragment velocity", " ", fi.velocity)
                 Pkfi = kill_probability(fi)
                 println("pki for fragment", " ", Pkfi )
                 println("number of fragments", " ", n_dispi)
                 Pki = 1-(1-Pkfi)^n_dispi
                 Pkhi = Pki*Phiti
                 println(" pki for projectile", " ", Pki)
                 println(" pkhi for projectile", " ", Pkhi)
                 Pkn=Pkn*(1-Pki)
                 println(" pkill" , " ", 1-Pkn)
                 Pkhn=Pkhn*(1-Pkhi)

                 println(" pkill if hit", " ", 1-Pkhn)
             else
                 miss = miss +1
                 println("target missed")

             end
#                Pnhitold = Pnhit
#                Ekinbold = Ekinb
#                velocityold = velocity
                 burst=burst+1
                 i = i+1
                 println("burst", " ", burst)
                 println("missed", " ", miss, " ", "times")
            end
            Pktot = (1-Pkhn)
            #Pktot = 1-Pkn
            #Ekinb = burst * 0.5*mass*(velocity)^2
            #burst, Pktot, Ekinb, mean(VEND)
            #@label End
            burst, Pktot, Ekinb


    end

end
