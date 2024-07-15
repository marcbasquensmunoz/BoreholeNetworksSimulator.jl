
abstract type TimeSuperpositionMethod end

# Compute the coefficient matrix of the ground model equations
function ground_model_coeffs!(M, ::TimeSuperpositionMethod, borefield) end
# Compute the independent vector of the ground model equations
function ground_model_b!(b, ::TimeSuperpositionMethod, borefield, step) end
# Compute auxiliary quantities needed for this particular model
function precompute_auxiliaries!(::TimeSuperpositionMethod, borefield, t) end
# Update auxiliary quantities after each time step
function update_auxiliaries!(::TimeSuperpositionMethod, X, borefield, step) end