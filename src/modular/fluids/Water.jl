
struct Water <: Fluid end

cpf(::Water) = 4182.
thermophysical_properties(::Water, Tref) = thermophysical_properties(Tref, "Water")
