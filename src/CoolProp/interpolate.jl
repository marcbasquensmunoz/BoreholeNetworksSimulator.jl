
function interpolate(T, property, value)
    i = 1
    while i < length(T) && T[i] < value
        i += 1
    end
    if i == length(T)
        return property[i]
    end
    return (property[i] + property[i+1])/2
end

function evaluate_thermophysical_properties(properties::ThermophysicalProperties, T)
    μ  = interpolate(properties.T, properties.μ,  T)
    ρ  = interpolate(properties.T, properties.ρ,  T)
    cp = interpolate(properties.T, properties.cp, T)
    k  = interpolate(properties.T, properties.k,  T)
    return (μ, ρ, cp, k)
end
