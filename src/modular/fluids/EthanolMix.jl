
"""
    EthanolMix <: Fluid (
        stored_properties::ThermophysicalProperties{Float64}
    )

Models a 20% ethanol and water mix.

To initialize, use the convenience method:
    function EthanolMix()
that will automatically compute `stored_properties`.
"""
struct EthanolMix <: Fluid 
    stored_properties::ThermophysicalProperties{Float64}
end
function EthanolMix()
    EthanolMix(load_properties("INCOMP::MEA-20%"))
end

cpf(::EthanolMix) = 4182.
thermophysical_properties(f::EthanolMix, T) = evaluate_thermophysical_properties(f.stored_properties, T)
