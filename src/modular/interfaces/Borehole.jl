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


"""
    BoreholeMock <: Borehole 

Mock for testing purposes.
"""
@with_kw struct BoreholeMock <: Borehole 
    H = 0.
    D = 0.
    rb = 0.
    Ci = 0.
    Co = 0.
    Cb = 0.
end
get_H(bh::BoreholeMock) = bh.H
get_D(bh::BoreholeMock) = bh.D
get_rb(bh::BoreholeMock) = bh.rb
uniform_Tb_coeffs(bh::BoreholeMock, λ, mass_flow, Tref, fluid) = bh.Ci, bh.Co, bh.Cb
