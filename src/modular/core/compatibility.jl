abstract type Compatibility end

struct Compatible <: Compatibility end
struct NotCompatible <: Compatibility
    message
end

check_compatibility(::Borefield, ::Constraint, ::Method) = Compatible()
