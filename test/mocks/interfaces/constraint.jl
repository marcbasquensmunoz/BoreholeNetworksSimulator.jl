
"""
    ConstraintMock <: Constraint 

Mock for testing purposes.
"""
@with_kw struct ConstraintMock <: Constraint 
    M = []
    b = []
end
BoreholeNetworksSimulator.constraints_coeffs!(M, c::ConstraintMock, operation, borefield) = c.M .= M
BoreholeNetworksSimulator.constraints_b!(b, c::ConstraintMock, operation, step) = c.b .= b
