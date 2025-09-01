abstract type Compatibility end

struct Compatible <: Compatibility end
struct NotCompatible <: Compatibility
    message
end

check_compatibility(::Medium, ::Constraint, ::TimeSuperpositionMethod) = Compatible()

function check_compatibility(medium::Medium, ::Constraint, ::OriginalNonHistoryMethod) 
    if medium isa FlowInPorousMedium 
        NotCompatible("The non-history method is not implemented with ground water flow yet") 
    else
        Compatible() 
    end
end