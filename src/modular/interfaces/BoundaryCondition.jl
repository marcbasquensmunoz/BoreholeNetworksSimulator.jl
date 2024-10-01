"""
    abstract type BoundaryCondition

Interface for boundary conditions.
"""
abstract type BoundaryCondition end



"""
    BoundaryConditionMock <: BoundaryCondition 

Mock for testing purposes.
"""
struct BoundaryConditionMock <: BoundaryCondition end
coefficients(::BoundaryConditionMock, setup::SegmentToSegment, params::Constants, dp, containers) = ones(length(containers.aux))
q_coef(::BoundaryConditionMock, m, method, sts, λ, i) = m.q_coef
