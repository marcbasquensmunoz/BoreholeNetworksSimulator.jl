
struct EthanolMix <: Fluid 
    #percentage = 0
end

cpf(::EthanolMix) = 4182.
thermophysical_properties(::EthanolMix, Tref) = thermophysical_properties(Tref, "INCOMP::MEA-20%")
