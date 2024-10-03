"""
    abstract type Borehole

Interface for boreholes.

Required functions:
- `get_H(::Borehole)`: Return the length of the borehole.
- `get_D(::Borehole)`: Return the burial depth of the borehole.
- `get_rb(::Borehole)`: Return the radius of the borehole.
- `uniform_Tb_coeffs(::Borehole, λ, mass_flow, Tref, fluid)`: Return the internal model coefficients
    for the resistance network between the pipes and the wall.
"""
abstract type Borehole end

@required Borehole begin
    get_H(::Borehole)
    get_D(::Borehole)
    get_rb(::Borehole)
    uniform_Tb_coeffs(::Borehole, λ, mass_flow, Tref, fluid)
end

get_n_segments(::Borehole) = 1
get_h(bh::Borehole) = get_H(bh)
