
function compute_response!(g, medium::FlowInPorousMedium, borefield::Borefield, boundary_condition::BoundaryCondition, t) 
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
