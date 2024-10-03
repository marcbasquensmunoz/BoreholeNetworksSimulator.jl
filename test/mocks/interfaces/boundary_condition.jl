
"""
    BoundaryConditionMock <: BoundaryCondition 

Mock for testing purposes.
"""
struct BoundaryConditionMock <: BoundaryCondition end

#BoreholeNetworksSimulator.coefficients(::BoundaryConditionMock, setup, params, dp, containers) = ones(length(containers.aux))
BoreholeNetworksSimulator.q_coef(::BoundaryConditionMock, m, method, sts, λ, i) = m.q_coef
