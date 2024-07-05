
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