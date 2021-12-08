function aimpoint(target::TargetRect)
    return (target.a/2,target.b/2)
end
"""
    adjustedFire(t, proj, weapon, aero,tank)

Adjsut the elevation and azimuth angle to hit the target within a precision of 1 mm. The function returns the adjusted
elevation and azimuth angles in degrees.
"""
function adjustedFire!(t::AbstractTarget, proj::AbstractPenetrator, weapon::Gun, aero::DataFrame,tank::Tank)
    precision=0.001
    missDist = 1e3
    QE = 0.0
    AZ = 0.0
    while missDist >precision
    QE,AZ=QEfinderMPMM(t, proj, weapon,aero)
    weapon.QE = QE#rad2deg(QE)
    weapon.AZ = AZ#rad2deg(AZ)
    tank.canon.θ = QE#rad2deg(QE)
    tank.turret.ξ = AZ#rad2deg(AZ)
    proj.position=muzzlePos(tank)
    proj.velocity=muzzleVel(tank)
    impactP= trajectoryMPMM(proj, t, weapon,aero)[1]
    missDist = euclidean(impactP,t.position)
    end 
    return QE,AZ
    
end 