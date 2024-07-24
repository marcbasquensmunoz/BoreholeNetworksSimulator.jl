
"""
    abstract type Medium

Common interface for mediums
"""
abstract type Medium end

# Get the thermal conductivity of the medium
function get_λ(::Medium) end
# Get the thermal diffusivity of the medium
function get_α(::Medium) end
# Compute the thermal responses between segments for all times in t 
function compute_response!(::Medium, borefield, coord_source, coord_eval, t) end