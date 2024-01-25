"""
precompute_auxiliaries!(model, borefield, t)               Compute auxiliary quantities needed for this particular model
update_auxiliaries!(model, X, borefield, step)             Update auxiliary quantities after each time step

ground_model_coeffs!(M, model, borefield)                  Compute the coefficient matrix of the ground model equations
ground_model_b!(b, model, borefield, step)                 Compute the independent vector of the ground model equations
"""
abstract type GroundModel end
