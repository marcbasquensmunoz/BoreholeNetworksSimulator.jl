
function load_properties(fluid_name)
    path = "$(@__DIR__)/../modular/fluids/data_$fluid_name"
    df = CSV.read(path, DataFrame)
    ThermophysicalProperties(T = Vector(df.T), μ = Vector(df.μ), ρ = Vector(df.ρ), cp = Vector(df.cp), k = Vector(df.k))
end