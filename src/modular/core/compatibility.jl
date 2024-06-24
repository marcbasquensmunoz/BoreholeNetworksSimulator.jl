abstract type Compatibility end

struct Compatible <: Compatibility end
struct NotCompatible <: Compatibility
    message
end

check_compatibility(borefield::Borefield, constraint::Constraint, method::Method) = Compatible()
