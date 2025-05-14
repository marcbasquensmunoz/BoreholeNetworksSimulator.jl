
global const fluid_names = Dict([
    :water => "Water",
    :ethanol20 => "EthanolMix20",
    :glycol30 => "GlycolMix30",
])

global const fluid_names_coolprop = Dict([
    :water => "Water",
    :ethanol20 => "INCOMP::MEA-20%",
    :glycol30 => "INCOMP::MEG-30%",
])

function load_properties(fluid_name)
    path = joinpath(dirname(@__DIR__), "modular", "fluids", "data_$fluid_name")
    df = CSV.read(path, DataFrame)
    ThermophysicalProperties{Float64}(T = Vector{Float64}(df.T), μ = Vector{Float64}(df.μ), ρ = Vector{Float64}(df.ρ), cp = Vector{Float64}(df.cp), k = Vector{Float64}(df.k))
end