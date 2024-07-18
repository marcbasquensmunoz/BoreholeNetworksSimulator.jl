abstract type Compatibility end

struct Compatible <: Compatibility end
struct NotCompatible <: Compatibility
    message
end

check_compatibility(::Borefield, ::Constraint, ::TimeSuperpositionMethod) = Compatible()

function check_compatibility(borefield::Borefield, ::Constraint, ::NonHistoryMethod) 
    if borefield.medium isa FlowInPorousMedium 
        NotCompatible("The non-history method is not implemented with ground water flow yet") 
    else
        Compatible() 
    end
end