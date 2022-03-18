function zoneVector(z::ZoneZdata, proj::AbstractPenetrator)
    
    
    velLower= [z.velocity[1]*cosd(z.angle[1]),z.velocity[1]*sind(z.angle[1]), 0.0]+[norm(proj.velocity),0.0,0.0]
    
    
    angleLower= acosd(clamp(dot([1.0,0.0,0.0],velLower)/(norm([1.0,0.0,0.0])*norm(velLower)), -1, 1))
    
    
    velUpper= [z.velocity[3]*cosd(z.angle[3]),z.velocity[3]*sind(z.angle[3]), 0.0]+[norm(proj.velocity),0.0,0.0]
    
    
    angleUpper=acosd(clamp(dot([1.0,0.0,0.0],velUpper)/(norm([1.0,0.0,0.0])*norm(velUpper)), -1, 1))
    
    return angleLower, angleUpper
    end 
    
    function inCone(proj::AbstractPenetrator, target::AbstractTarget,bounds::Tuple{Float64,Float64})
        r = target.position - proj.position
        
     angle=acosd(clamp(dot(proj.velocity,r)/(norm(proj.velocity)*norm(r)), -1, 1))
    
    
    if angle > bounds[1] && angle < bounds[2]
        in= true  
    else
        in=false
    end  
    return in
    end 
    
    function hit(zdata::ZoneList,proj::AbstractPenetrator,target::AbstractTarget)
        inside = []
        for z in zdata.list
            bounds = zoneVector(z,proj)
            push!(inside,inCone(proj,target,bounds))
    
        end 
        return inside
    end 
#end
