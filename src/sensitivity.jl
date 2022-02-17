function impactAbm(abm::AbstractPenetrator, tVulnerable::CSGGenerator{Float64})

    fragments = frag(abm, zones,shapes)
 
    frs,rays = raytracing(fragments)
 
     
    hit,shots = shotlines(tVulnerable,rays)
    nfrs = length(frs)
 

 return hit,frs
 
end 