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
coefficients_sts(::BoundaryConditionMock, sts::SegmentToSegment, params::Constants, dp) = 1.