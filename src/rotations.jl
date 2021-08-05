

function targetPos(target::Target, tourelle::Tourelle, gun::Canon)
rdoelprimprim = [target.ρ; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tourelle.Φ), -deg2rad(tourelle.ψ), deg2rad(tourelle.χ), :XZY)
RLos = angle_to_dcm(0, -deg2rad(gun.θₜ), deg2rad(gun.ξₜ), :XZY)
rdoelprim = Rhull*RLos*rdoelprimprim
rdoel = rdoelprim + [0; tourelle.rWY; 0]
return rdoel
end

function muzzlePos(gun::Canon, tourelle::Tourelle )
rprimprim = [gun.lw; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tourelle.Φ), -deg2rad(tourelle.ψ), deg2rad(tourelle.χ), :XZY)
Rloop = angle_to_dcm(0, -deg2rad(tourelle.θ), 0, :XZY)
Rtoren = angle_to_dcm(0, 0, deg2rad(tourelle.ξ), :XZY)
rprim = Rhull*Rtoren*Rloop*rprimprim
r = rprim +[0; tourelle.rWY; 0]
return r
end

function muzzleVel(gun::Canon,tourelle::Tourelle )
uprimprim = [gun.u₀; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tourelle.Φ), -deg2rad(tourelle.ψ), deg2rad(tourelle.χ), :XZY)
Rloop = angle_to_dcm(0, -deg2rad(tourelle.θ), 0, :XZY)
Rtoren = angle_to_dcm(0, 0, deg2rad(tourelle.ξ), :XZY)
uprim = Rhull*Rtoren*Rloop*uprimprim
u =uprim
return u
end

function wind(tourelle::Tourelle,canon::Canon, w::Wind)
    w = [-w.cross*sind(tourelle.ξ+canon.ξₜ)-w.range*cosd(tourelle.ξ+canon.ξₜ);0.0;w.cross*cosd(tourelle.ξ+canon.ξₜ)-w.range*sind(tourelle.ξ+canon.ξₜ)]
end
