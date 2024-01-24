using Parameters
using GeometryTypes
using BoreholeResponseFunctions

@with_kw struct GroundWaterMedium <: Medium
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

λ(bfp::GroundWaterMedium) = bfp.λ

function response(medium::GroundWaterMedium, borefield::Borefield, coord_source, coord_eval, t) 
    p =  GeometryTypes.Point3{Float64}.(coord_source)
    tp = GeometryTypes.Point3{Float64}.(coord_eval)  

    # Rotation of points in new coordinate system where Darcy velocity is parallel to x axis
    p_rot  = rotation_z(p, -medium.θ) 
    tp_rot = rotation_z(tp, -medium.θ) 

    distances = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp_rot) 
    d = [d[1] == 0. && d[2] == 0. ?  (0., get_rb(borefield, i), d[3], d[4]) : d for (i, d) in enumerate(distances)]

    g = [1/(2π*medium.λs)*mfls_adiabatic_surface(tt, medium.α, coord[1:3]..., medium.vt, get_h(borefield, i), coord[4]; atol =1e-9) for (i, coord) in enumerate(d), tt in t]
    return g
end