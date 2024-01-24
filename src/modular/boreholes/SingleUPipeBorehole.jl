using Parameters

@with_kw struct SingleUPipeBorehole <: Borehole
    λg = 1.5                            # grout conductivity
    Cg = 2000. * 1550.                  # grout capacity
    αg = λg/Cg	                        # grout thermal diffusivity

    rp = 0.02                           # equivalent pipe radius
    λp = 0.42                           # p - pipe material
    dpw = 0.0023                        # pipe thickness
    rpo = rp + dpw                      # equivalent pipe radius
    hp = 725.                           # heat transfer coefficient fluid to pipe ?
    pipe_position = Point2.([[0.03,0.0], [-0.03,.0]]) 
        
    rb = 0.115/2                        # borehole radius
    H                                   # length of the borehole
    D                                   # burial depth of the borehole
end

get_H(bh::SingleUPipeBorehole) = bh.H
get_D(bh::SingleUPipeBorehole) = bh.D
get_h(bh::SingleUPipeBorehole) = bh.H
get_rb(bh::SingleUPipeBorehole) = bh.rb

function resistance_network(borehole::SingleUPipeBorehole, λs, mass_flow)
    x = [p.data[1] for p in borehole.pipe_position]
    y = [p.data[2] for p in borehole.pipe_position]
    
    @unpack λg, λp, rb, rp, rpo, dpw, hp = borehole

    # Make hp as a function of mass flow
    Rp = 1/(2*pi*λp)*log(rp/(rp-dpw)) + 1/(2*pi*rp*hp)   # pipe resistance
    
    N = length(borehole.pipe_position)
    R = zeros(N, N)

    for i in 1:N
        for j in 1:N
            if i == j
                R[i,j] = 1/(2*pi*λg) * ( log(rb/rpo) - (λg - λs)/(λg + λs) * log(1 - (x[j]^2 + y[j]^2) / rb^2) ) + Rp
            else
                dij = sqrt( (1 - (x[i]^2+y[i]^2) / rb^2) * (1 - (x[j]^2+y[j]^2) / rb^2) + ( (x[i] - x[j])^2 + (y[i] - y[j])^2) / rb^2 )
                R[i,j] = - 1/(2*pi*λg) * (log(( (x[i] - x[j])^2 + (y[i] - y[j])^2) / rb^2 ) + (λg - λs)/(λg + λs) * log(dij))
            end
        end
    end
    return R
end
