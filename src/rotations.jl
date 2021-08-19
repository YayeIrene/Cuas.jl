

function targetPos(target::AbstractTarget, tank::Tank)
rdoelprimprim = [target.ρ; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tank.hull.Φ), -deg2rad(tank.hull.ψ), deg2rad(tank.hull.χ), :XZY)
RLos = angle_to_dcm(0, -deg2rad(tank.sight.θ), deg2rad(tank.sight.ξ), :XZY)
rdoelprim = Rhull*RLos*rdoelprimprim
rdoel = rdoelprim + [0; tank.rWY; 0]
return rdoel
end

function muzzlePos(tank::Tank)
rprimprim = [tank.canon.lw; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tank.hull.Φ), -deg2rad(tank.hull.ψ), deg2rad(tank.hull.χ), :XZY)
Rloop = angle_to_dcm(0, -deg2rad(tank.canon.θ), 0, :XZY)
Rtoren = angle_to_dcm(0, 0, deg2rad(tank.turret.ξ), :XZY)
rprim = Rhull*Rtoren*Rloop*rprimprim
r = rprim +[0; tank.rWY; 0]
return r
end

function muzzleVel(tank::Tank)
uprimprim = [tank.canon.u₀; 0; 0]
Rhull = angle_to_dcm(-deg2rad(tank.hull.Φ), -deg2rad(tank.hull.ψ), deg2rad(tank.hull.χ), :XZY)
Rloop = angle_to_dcm(0, -deg2rad(tank.canon.θ), 0, :XZY)
Rtoren = angle_to_dcm(0, 0, deg2rad(tank.turret.ξ), :XZY)
uprim = Rhull*Rtoren*Rloop*uprimprim
u =uprim
return u
end

function wind(tank::Tank, w::Wind)
    w = [-w.cross*sind(tank.turret.ξ+tank.hull.χ)-w.range*cosd(tank.turret.ξ+tank.hull.χ);0.0;w.cross*cosd(tank.turret.ξ+tank.hull.χ)-w.range*sind(tank.turret.ξ+tank.hull.χ)]
end
