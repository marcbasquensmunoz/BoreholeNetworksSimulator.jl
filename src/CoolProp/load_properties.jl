
function load_properties(fluid_name)
    path = joinpath(dirname(@__DIR__), "modular", "fluids", "data_$fluid_name")
    df = CSV.read(path, DataFrame)
    ThermophysicalProperties{Float64}(T = Vector{Float64}(df.T), μ = Vector{Float64}(df.μ), ρ = Vector{Float64}(df.ρ), cp = Vector{Float64}(df.cp), k = Vector{Float64}(df.k))
end