using CoolProp
using CSV
using DataFrames

fluid_name = "Water"
T_lims = (0, 100)


pressure = 101325
T0 = 273.153

T_delta = 0.1
T_sample = T_lims[1]:T_delta:T_lims[2]

μ = zeros(length(T_sample))
ρ = zeros(length(T_sample))
cp = zeros(length(T_sample))
k = zeros(length(T_sample))

for (i, T) in enumerate(T_sample)
    T_ref = T + T0
    μ[i]  = PropsSI("viscosity",    "T", T_ref, "P", pressure, fluid_name)
    ρ[i]  = PropsSI("D",            "T", T_ref, "P", pressure, fluid_name)
    cp[i] = PropsSI("C",            "T", T_ref, "P", pressure, fluid_name)
    k[i]  = PropsSI("conductivity", "T", T_ref, "P", pressure, fluid_name)	
end

df = DataFrame(T = T_sample, μ = μ, ρ = ρ, cp = cp, k = k)
CSV.write("$(@__DIR__)/../modular/fluids/data_$fluid_name", df)