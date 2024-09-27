
"""
    Water <: Fluid (
        stored_properties::ThermophysicalProperties{Float64}
    )

Models water.

To initialize, use the convenience method:
    function Water()
that will automatically compute `stored_properties`.
"""
struct Water <: Fluid 
    stored_properties::ThermophysicalProperties{Float64}
end
function Water()
    Water(load_properties("Water"))
end

cpf(::Water) = 4182.
thermophysical_properties(f::Water, T) = evaluate_thermophysical_properties(f.stored_properties, T)
