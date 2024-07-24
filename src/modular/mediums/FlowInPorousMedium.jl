"""
    FlowInPorousMedium{T <: Real} <: Medium @deftype T

Model a porous ground with a water flow.

# Arguments
- `λw = 0.6`: water thermal conductivity
- `λs = 2.`: ground thermal conductivity
- `Cw = 4.18*1e6`:water thermal capacity
- `Cs = 1.7*1e6`: ground thermal capacity
- `θ = 0.`: angle of Darcy velocity
- `Φ = 0.2 `: porosity
- `λ = λs * (1-Φ) + λw*Φ`: porous medium conductivity
- `C = Cs * (1-Φ) + Cw*Φ`: porous medium capacity
- `α = λ/C`: porous medium thermal diffusivity
- `ux_in_meterperday = 1e-2`: groundwater speed along the flow coordinate
- `ux = ux_in_meterperday/(3600*24)`: groundwater speed in m/s
- `vt = ux * Cw/C`: porous medium darcy velocity
- `T0 = 0.`: initial ground temperature
"""
@with_kw struct FlowInPorousMedium{T <: Real} <: Medium @deftype T
    λw = 0.6                            # water thermal conductivity
    λs = 2.                             # ground thermal conductivity
    Cw = 4.18*1e6                       # water thermal capacity
    Cs = 1.7*1e6                        # ground thermal capacity
    θ = 0.                              # angle of Darcy velocity
    Φ = 0.2                             # porosity
    λ = λs * (1-Φ) + λw*Φ               # porous medium conductivity
    C = Cs * (1-Φ) + Cw*Φ               # porous medium capacity
    α = λ/C                             # porous medium thermal diffusivity

    ux_in_meterperday = 1e-2            # groundwater speed along the flow coordinate
    ux = ux_in_meterperday/(3600*24)    # groundwater speed in m/s
    vt = ux * Cw/C                      # porous medium darcy velocity

    T0 = 0.
end

get_λ(medium::FlowInPorousMedium) = medium.λ
get_α(medium::FlowInPorousMedium) = medium.α
get_T0(medium::FlowInPorousMedium) = medium.T0
