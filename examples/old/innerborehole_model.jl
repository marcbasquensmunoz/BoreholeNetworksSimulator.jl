# DEFAULTS for borhole parameters
@with_kw struct BoreholePara{R<:Real} @deftype R
    
    rb = 0.115/2 

    #Grout    
    λg = 1.5     
    Cg = 2000. * 1550.          
    αg = λg/Cg	

    # Ground
    λs = 3.
    Cs = 2500. * 750.           
    αs = λs/Cs

    #Pipe Fluid
    # rp = 0.028*sqrt(2)          #equivalent pipe radius
    # rp = 0.02*sqrt(2)           #equivalent pipe radius
    rp = 0.02                   #equivalent pipe radius
    λp = 0.42;                  # p - pipe material
    dpw = 0.0023;               # pipe thickness
    hp = 725.;                  # heat transfer coefficient water to pipe ?
    rpo = rp + dpw                   #equivalent pipe radius
    
    # Cw = pi*rp^2 *4182*1000.    # Water capacity
    
    Rp = 1/(2*pi*λp)*log(rp/(rp-dpw)) + 1/(2*pi*rp*hp)   # pipe resistance


    #Resistances
    # Rgrout = 1/(2*pi*λb)*log(rb/rp)                  # grout resistance
    # Rb = (Rgrout + Rp)     # in this example the borehole resistance is considered         
end


"""
Use the line source approximation to calculate the resistance matrix
on a borehole section.
"""
function resistance_network(params::BoreholePara, positions)
    x = [p.data[1] for p in positions]
    y = [p.data[2] for p in positions]
    
    @unpack λg,λs,rb,rpo,Rp = params

    N = length(positions)
    R = zeros(N,N)

    for i = 1:N
        for j = 1:N
            if i==j
                R[i,j] = 1/(2pi*λg)* ( log(rb/rpo) - (λg - λs)/(λg + λs) * log(1 - (x[j]^2 + y[j]^2) / rb^2) ) + Rp
            else
                dij = sqrt( (1 - (x[i]^2+y[i]^2) / rb^2 ) * (1 - (x[j]^2+y[j]^2) / rb^2 )  + 
                            ( (x[i] - x[j])^2 + (y[i] - y[j])^2) / rb^2  )
                R[i,j] =  -1/(2pi*λg) * (log(( (x[i] - x[j])^2 + (y[i] - y[j])^2) / rb^2  ) + (λg - λs)/(λg + λs) * log(dij))
            end

        end
    end
    return R
end


"""
Return RΔ1,RΔ2,RΔ12 for a given borehole configuration using line source model 
(first order multipole approximation)
"""
function deltacircuit(params::BoreholePara, positions)
        
    if length(positions) != 2
         throw("This function is designed only for single U-pipes")
         return nothing
    end

    R = resistance_network(params, positions) 
    
    RΔ1  = (R[1,1]*R[2,2] - R[1,2]^2)/(R[2,2] - R[1,2])
    RΔ2  = (R[1,1]*R[2,2] - R[1,2]^2)/(R[1,1] - R[1,2])  
    RΔ12 = (R[1,1]*R[2,2] - R[1,2]^2)/(R[1,2]) 

    return RΔ1,RΔ2,RΔ12
end

"""
Return the effective borehole resistance Rbstart for uniform borehole temperature 
boundary condition 
"""
function effective_borehole_resistance(params::BoreholePara, positions, H, Cf, Vf)
    
    RΔ1,RΔ2,RΔ12 =  deltacircuit(params, positions)
    Rb = (RΔ1 * RΔ2)/(RΔ1 + RΔ2)
    η = H/Cf*Vf * 1/2Rb * sqrt(1 + 4*Rb/RΔ12)
    Rbstar  = Rb * η *coth(η)
    return Rbstar
end


"""
Matrix A according Cimmino 2015
"""
function coefficient_matrix(R, Cf, Vf)

    A = inv(R)
    N = size(R)[1]
    np = div(N,2)
    
    for j = 1:N
        for i = 1:N
            direction = i <= np ? -1 : 1
            A[i,j] *= direction / (Cf*Vf)
        end
    end
    
    return A
end


"""
(Cimmino 2016)
"""
function uniformTb_koeff(A, H)        
    N = size(A)[1]
    np = div(N, 2)

    EH = exp(A*H)
    @views EoutH = EH[1:np , np+1:2np] - EH[np+1:2np,np+1:2np]
    @views EinH  = EH[np+1:2np , 1:np] - EH[1:np,1:np]      
    
    return EinH, -EoutH, (EoutH - EinH)*ones(1,np)
end
