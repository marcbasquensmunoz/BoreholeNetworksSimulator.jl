
"""
    GlycolMix <: Fluid (
        stored_properties::ThermophysicalProperties{Float64}
    )

Models a 30% glycol and water mix.

To initialize, use the convenience method:
    function GlycolMix()
that will automatically compute `stored_properties`.
"""
struct GlycolMix <: Fluid 
    stored_properties::ThermophysicalProperties{Float64}
end
function GlycolMix()
    GlycolMix(load_properties(fluid_names[:glycol30]))
end

cpf(::GlycolMix) = 3710.
thermophysical_properties(f::GlycolMix, T) = evaluate_thermophysical_properties(f.stored_properties, T)
