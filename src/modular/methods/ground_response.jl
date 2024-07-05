
function compute_response!(g, medium::GroundMedium, borefield::Borefield, t) 
    @unpack λ, α = medium
    Ns = segment_amount(borefield)

    coords = [segment_coordinates(borefield, i) for i in 1:Ns]
    
    for (k, tt) in enumerate(t)
        for j in 1:Ns
            for i in 1:Ns
                x1, y1, zref1, _  = coords[i]
                x2, y2, zref2, H2 = coords[j]
                zeval2 = zref2 + H2/2
                @show coords[i]
                @show coords[j]
                @show zeval2


                x = x2 - x1
                y = abs(y2 - y1) 
                D = zref1
                z = zeval2
                if x == 0. && y == 0.
                    y = get_rb(borefield, i)
                end
                h = get_h(borefield, i)
                #g[i, j, k] = fls_dirichlet(tt, x, y, z, h, D, α, λ) 
                g[i, j, k] = fls(tt, x, y, z, h, D, α, λ) 
            end
        end
    end
end

point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)    
fls_step_response(t, σ, z, a, b, α, kg; atol=1e-8) = quadgk(ζ -> point_step_response(t, sqrt(σ^2 + (z-ζ)^2), α, kg), a, b, atol=atol)[1]

function fls_adiabatic_surface(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    I1 = fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
    I2 = fls_step_response(t, σ, z, -D-H, -D, α, kg, atol=atol) 
    return (I1 + I2)
end

function fls_dirichlet(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    I1 = fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
    I2 = fls_step_response(t, σ, z, -D-H, -D, α, kg, atol=atol) 
    return (I1 - I2)
end

function fls(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    return fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
end
