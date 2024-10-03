
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
BoreholeNetworksSimulator.get_H(bh::BoreholeMock) = bh.H
BoreholeNetworksSimulator.get_D(bh::BoreholeMock) = bh.D
BoreholeNetworksSimulator.get_rb(bh::BoreholeMock) = bh.rb
BoreholeNetworksSimulator.uniform_Tb_coeffs(bh::BoreholeMock, λ, mass_flow, Tref, fluid) = bh.Ci, bh.Co, bh.Cb
