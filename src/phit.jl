function dispersionGrid(ϵ::RectErrorB, pmi::Vector{Float64};N = 100000, xmin=-0.5, ymin=-0.5,xmax=0.5,ymax=0.5,dim=10,sizex=2.3,sizey=2.3)
    d2 = MvNormal([ϵ.μ_x,ϵ.μ_y],[(ϵ.σ_x)^2 0.0 ;0.0 (ϵ.σ_y)^2])
    #N=100000
    x2 = rand(d2,N)
    p=zeros(N)
    #count(i->(4.2<=i<=6), [2,3,4.1,5.2,6])

    #X = xmin+pmi[1]:0.1:xmax+pmi[1]
    #Y = ymin+pmi[2]:0.1:ymax+pmi[2]
    #g = Iterators.product(X,Y)

    #collect.(p)
    #grid = vec(collect.(g))
    #dim = 10
    g = CartesianGrid((dim,dim),(-sizex/2+pmi[1],-sizey/2+pmi[2]),(sizex/dim,sizey/dim))
    p = collect(elements(g))
    impact = Array{Point2}(undef,N)
    for i=1:N
        impact[i] = Point(x2[1,i],x2[2,i])
    end 
    hit= zeros(length(p))
    cg = []
    for i=1:length(p)
        for j=1:N
        if (impact[j]∈ p[i])
            hit[i]+=1
            #push!(cg,impact[j])
        end
        end
    end 

    phit=hit./N
    

    grid = reshape(phit,dim,dim)

    for cell in p 
        centroid(cell)
    end

    
end 
