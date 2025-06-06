using CoolProp
using CSV
using DataFrames

function generate_data(fluid, T_lims)
    pressure = 101325
    T0 = 273.153

    T_delta = 0.1
    T_sample = T_lims[1]:T_delta:T_lims[2]

    μ = zeros(length(T_sample))
    ρ = zeros(length(T_sample))
    cp = zeros(length(T_sample))
    k = zeros(length(T_sample))

    fluid_name = fluid_names_coolprop[fluid]
    for (i, T) in enumerate(T_sample)
        T_ref = T + T0
        μ[i]  = PropsSI("viscosity",    "T", T_ref, "P", pressure, fluid_name)
        ρ[i]  = PropsSI("D",            "T", T_ref, "P", pressure, fluid_name)
        cp[i] = PropsSI("C",            "T", T_ref, "P", pressure, fluid_name)
        k[i]  = PropsSI("conductivity", "T", T_ref, "P", pressure, fluid_name)	
    end

    df = DataFrame(T = T_sample, μ = μ, ρ = ρ, cp = cp, k = k)
    path = joinpath(dirname(@__DIR__), "modular", "fluids", "data_$(fluid_names[fluid])")
    CSV.write(path, df)
end

generate_data(:water, (0, 100))
generate_data(:ethanol20, (-11.118, 40))
generate_data(:glycol30, (-11.118, 40))