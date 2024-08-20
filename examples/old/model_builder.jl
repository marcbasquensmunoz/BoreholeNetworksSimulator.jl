# BUILD OF THE MODEL 
#example for 2 sources
# 8 unknowns Tfin1, Tfout1, Tfin2, Tfout2, Tb1, Tb2, qb1, qb2
# If we impose the temperature Tin this problem can be solved in the Laplace domain

# INJECTION 
#  _                                                                          _   _         _       _                           _
# |  kin    kout    0      0      kb1    0       0         0       0       0   |  |  Tf_in1  |     |  0                          |  1. Borehole Model  
# |  0      0       kin    kout   0      kb2     0         0       0       0   |  |  Tf_out1 |     |  0                          |  
# |  1      0       0      0      0      0       0         0       0       0   |  |  Tf_in2  |     |  Tin1(t)                    |  2. Constrain on input temperature
# |  0     -1       1      0      0      0       0         0       0       0   |  |  Tf_out2 |  =  |  0                          |
# |  0      0       0      0   -1*2πkg   0      g11      g12       g13     g14 |  |  Tb1     |     |  -T0 - effect_of_past_loads |  3. Ground Model
# |  0      0       0      0   -1*2πkg   0      g21      g22       g23     g24 |  |  Tb2     |     |  -T0 - effect_of_past_loads |
# |  0      0       0      0      0   -1*2πkg   g31      g32       g33     g34 |  |  Δqb1'   |     |  -T0 - effect_of_past_loads |  
# |  0      0       0      0      0   -1*2πkg   g41      g32       g43     g44 |  |  Δqb2'   |     |  -T0 - effect_of_past_loads |
# | mcp   -mcp      0      0      0      0      -h        -h       0       0   |  |  Δqb3'   |     | +H*sum(Δqb1'past)           |  4. Heat Balance
# |_ 0      0      mcp   -mcp     0      0       0         0       -h     -h  _|  |_ Δqb4'  _|     |_+H*sum(Δqb2'past)          _|



"""
build_matrix for uniform borehole temperature but allowing the use of multiple segments per borehole
M  - matrix 
Nb - number of borehole 
Ns - segments
"""
function build_matrix!(M,Nb,Ns,
    k_in, k_out, k_b, 
    G,
    mf, cpf, bh_map, h,             
    branches
    )

    # BOREHOLE MODEL
    for i = 1:Nb    
        M[i,i*2-1]     =  k_in[1]
        M[i,i*2]       =  k_out[1]
        M[i,Nb*2 + i]   =  k_b[1]    
    end

    # CONSTRAIN ON INPUT TEMPERATURE 
    #We always assign the input temperature in each borehole
    for i = Nb+1:2Nb
        M[i,(i-Nb)*2-1] = 1
    end

    #branches constraints
    for br in branches      
        for i=2:length(br)
            M[Nb+br[i],(br[i-1])*2] = -1 
        end
    end

    #ground model
    M[2Nb+1:2Nb+Ns,3Nb+1:3Nb+Ns] = G
    for (j,i) in enumerate(2Nb+1:2Nb+Ns)
        M[i, 2Nb + bh_map[j] ] = -1  
    end

    # HEAT BALANCE
    for i=1:Nb
        M[2Nb+Ns+i , i*2-1:i*2] = [+mf*cpf -mf*cpf]
    end

    for (j,i) in enumerate(2Nb+Ns+1:3Nb+Ns)
        for k=1:Ns
            if bh_map[k] == j 
                M[i,k+3Nb] = -h
            # else
            #     M[i,k+3Nb] = 0.
            end
        end
    end   
end


function build_giventerm!(b,Nb,Ns,Tin,T0,branches)

    #branches constraints
    for br in branches              
        b[Nb+br[1]] = Tin      
        for i=2:length(br)
            b[Nb+br[i]] = 0.
        end
    end

    for i = 1:Ns
    b[2Nb+i] = -T0
    end
end


"""
Time domain scheme
update the given term b given nth_step, and previous loads.
"""

function update_b!(b, nth_step, Nb, Ns,
                    Tfin_constraint,branches,
                    qprime, g, T0,
                    Δqbcurrentsum, h, 
                    )

    b[1:end] .= 0.

    #branches constraints
    for br in branches              
        b[Nb+br[1]] = Tfin_constraint[nth_step]
    end

    for (idx,i) in enumerate(2Nb+1:2Nb+Ns)
        b[i] = -T0 #0.  there was a mistake here
        for j=1:Ns
            for k=2:nth_step
                b[i] +=  - qprime[nth_step - k + 1, j] * g[j,idx, k]
            end
        end
    end

    for i=1:Nb
        b[2Nb+Ns+i] = Δqbcurrentsum[i]*h
    end  

    return nothing 
end



"""
one step of solution
"""
function solve_full_convolution_step!(X,A,b, nth_step, Nb, Ns,
                            Tfin_constraint,branches,
                            qprime, g, T0,
                            Δqbcurrentsum, h                  
                            )
  
    update_b!(b, nth_step, Nb, Ns,
             Tfin_constraint, branches,
             qprime, g, T0,
             Δqbcurrentsum, h, 
            )
    x = A\b
    X[nth_step,:] = x

    qprime[nth_step,:] = x[3Nb+1:end]
    Δqbcurrentsum .+= x[3Nb+1:end]
end
