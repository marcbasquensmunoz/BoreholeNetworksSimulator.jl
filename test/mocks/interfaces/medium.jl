
"""
    MediumMock <: Medium

Mock for testing purposes.
"""
@with_kw struct MediumMock <: Medium 
    λ = 0.
    α = 0.
    T0 = 0.
    g = 0.
    q_coef = 0.
end
BoreholeNetworksSimulator.get_λ(m::MediumMock) = m.λ
BoreholeNetworksSimulator.get_α(m::MediumMock) = m.α
BoreholeNetworksSimulator.get_T0(m::MediumMock) = m.T0
BoreholeNetworksSimulator.compute_response!(g, m::MediumMock, options) = g .= m.g
BoreholeNetworksSimulator.constant_integral(::MediumMock, method, setup, λ, i) = 0.