
"""
    BorefieldMock <: Borefield 

Mock for testing purposes.
"""
@with_kw struct BorefieldMock <: Borefield 
    Nb = 1
    D = []
    H = []
    rb = []
    coordinates = []
    M = []
    b = []
end
BoreholeNetworksSimulator.n_boreholes(bf::BorefieldMock) = bf.Nb
BoreholeNetworksSimulator.get_D(bf::BorefieldMock, i) = bf.D[i]
BoreholeNetworksSimulator.get_H(bf::BorefieldMock, i) = bf.H[i]
BoreholeNetworksSimulator.get_rb(bf::BorefieldMock, i) = bf.rb[i]
BoreholeNetworksSimulator.segment_coordinates(bf::BorefieldMock, segment) = bf.coordinates[segment][1], bf.coordinates[segment][2], bf.D[segment], bf.H[segment]
BoreholeNetworksSimulator.internal_model_coeffs!(M, bf::BorefieldMock, medium, operation, T_fluid, fluid) = bf.M .= M
BoreholeNetworksSimulator.internal_model_b!(b, bf::BorefieldMock) = bf.b .= b
BoreholeNetworksSimulator.T_past_influence(bf::BorefieldMock) = zeros(bf.Nb)
