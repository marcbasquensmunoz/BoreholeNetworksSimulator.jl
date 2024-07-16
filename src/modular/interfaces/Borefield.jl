
abstract type Borefield end

# Get the amount of boreholes in the borefield     
function n_boreholes(::Borefield) end
# Get the total amount of segments in the borefield 
function n_segments(::Borefield) end
# Get the length of borehole i 
function get_H(::Borefield, i) end
# Get the segment length of borehole i
function get_h(::Borefield, i) end
# Get the radius of borehole i
function get_rb(::Borefield, i) end
# Get a vector containing to which borehole each segment belongs
function segment_map(::Borefield) end
# Get a vector with the coordinates of each segment  
function segment_coordinates(::Borefield) end
# Compute the coefficient matrix of the internal model equations
function internal_model_coeffs!(M, ::Borefield, operation, cpf) end
# Compute the independent vector of the internal model equations
function internal_model_b!(b, ::Borefield) end