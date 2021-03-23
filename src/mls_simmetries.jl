

struct GroundWaterFlow end
"""
return a vector of mutual distances pairs [(Δri,Δxi) for i = 1:N] 
"""
function evaluate_relevant_distances(::GroundWaterFlow, Ps1::Array{T},Ps2::Array{T}) where {T <: GeometryTypes.Point2}
    
    c1 = hcat(collect(Array(p) for p in Ps1)...)
    c2 = hcat(collect(Array(p) for p in Ps2)...)
    
    n,m = size(c1)[2], size(c2)[2]
    Δx = zeros(eltype(c1), (n,m))
    Δy = zeros(eltype(c2), (n,m))

    for j =1:m
        for k=1:n
            Δx[k,j] =  -c1[1,k] + c2[1,j]             
        end
    end

    for j =1:m
        for k=1:n
            Δy[k,j] =  abs(c1[2,k] - c2[2,j])
        end
    end 

    map(Δx,Δy) do x,y
        (x,y)
    end 
end


function evaluate_relevant_distances(::GroundWaterFlow, Ps1::Array{T},Ps2::Array{T}) where {T <: GeometryTypes.Point3}
    
    c1 = hcat(collect(Array(p) for p in Ps1)...)
    c2 = hcat(collect(Array(p) for p in Ps2)...)
    
    n,m = size(c1)[2], size(c2)[2]
    Δx = zeros(eltype(c1), (n,m))
    Δy = zeros(eltype(c1), (n,m))
    z = zeros(eltype(c1), (n,m))    
    D = zeros(eltype(c1), (n,m))    

    for j =1:m
        for k=1:n
            Δx[k,j] =  -c1[1,k] + c2[1,j]             
        end
    end

    for j =1:m
        for k=1:n
            Δy[k,j] =  abs(c1[2,k] - c2[2,j])
        end
    end 

    for j =1:m
        for k=1:n
            z[k,j] =  c2[3,j] 
            D[k,j] =  c1[3,k] 
        end
    end 

    map(Δx,Δy,z,D) do Δx,Δy,z,D
        (Δx,Δy,z,D)
    end 
end

"""
map unique pairs given a vector d of distance 
"""
function map_unique_pairs(d)
    du = unique(d)
    Amap   = fill(0,size(d));
    
    for i=1:length(d)
        Amap[i] =  findfirst(isequal(d[i]),du)
    end

    return du,Amap
end



# #example
# sources =  [(0.,1.),(3.,1.), (7.,1.), (0.,-1.), (3.,-1.), (7.,-1.), (0.,-4), (3.,-4.), (7.,-4.)]
# p = Point2{Float64}.(sources)  
# d = evaluate_relevant_distances(GroundWaterFlow(), p, p)
# du,m  = map_unique_pairs(d)   
# du[1] = (0.05,0.)

# vt = 1e-8
# α  = 1e-6
# t  = 1e8

# responses  = [mils(t, α, coord... ,  vt) for coord in du]

# solution =  responses[m]