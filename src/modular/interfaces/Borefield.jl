"""
borehole_amount                                                      Get the amount of boreholes in the borefield                
segment_amount                                                       Get the total amount of segments in the borefield 
get_H(i)                                                             Get the length of borehole i 
get_h(i)                                                             Get the segment length of borehole i
get_rb(i)                                                            Get the radius of borehole i
segment_map                                                          Get a vector containing to which borehole each segment belongs
segment_coordinates                                                  Get a vector with the coordinates of each segment    

internal_model_coeffs!(M, borefield, operation, cpf)                 Compute the coefficient matrix of the internal model equations
internal_model_b!(b, borefield)                                      Compute the independent vector of the internal model equations
"""
abstract type Borefield end

"""
get_λ                                                                Get the thermal conductivity of the medium
get_α                                                                Get the thermal diffusivity of the medium
compute_response!(medium, borefield, coord_source, coord_eval, t)    Compute the thermal responses between segments for all times in t                                 
"""
abstract type Medium end
