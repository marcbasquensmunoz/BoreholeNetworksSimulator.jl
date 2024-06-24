
@with_kw struct GroundWaterMedium{T <: Real} <: Medium @deftype T
    ux_in_meterperday = 1e-2            # groundwater speed along the flow coordinate
    ux = ux_in_meterperday/(3600*24)    # groundwater speed in m/s
    λw = 0.6                            # water thermal conductivity
    λs = 2.                             # ground thermal conductivity
    Cw = 4.18*1e6                       # water thermal capacity
    Cs = 1.7*1e6                        # ground thermal capacity
    θ = 0.                              # angle of Darcy velocity
    Φ = 0.2                             # porosity
    λ = λs *(1-Φ) + λw*Φ                # porous medium conductivity
    C = Cs *(1-Φ) + Cw*Φ                # porous medium capacity
    α = λ/C                             # porus medium thermal diffusivity
    vt = ux * Cw/C                      # porous medium darcy velocity
end

get_λ(medium::GroundWaterMedium) = medium.λ
get_α(medium::GroundWaterMedium) = medium.α

function compute_response!(g, medium::GroundWaterMedium, borefield::Borefield, t) 
    @unpack θ, λ, vt, α = medium
    Ns = segment_amount(borefield)

    coords = [segment_coordinates(borefield, i) for i in 1:Ns]
    
    for (k, tt) in enumerate(t)
        for j in 1:Ns
            for i in 1:Ns
            
                x1, y1, zref1, _  = coords[i]
                x2, y2, zref2, H2 = coords[j]
                zeval2 = zref2 + H2/2

                x = x2 - x1
                y = abs(y2 - y1) 
                D = zref1
                z = zeval2
                if x == 0. && y == 0.
                    y = get_rb(borefield, i)
                end
                h = get_h(borefield, i)
                g[i, j, k] = 1/(2π*λ)*mfls_adiabatic_surface(tt, α, x, y, z, vt, h, D; atol =1e-9) 
            end
        end
    end
end
