
function load_cache!(;containers::SimulationContainers, parameters::SimulationParameters, cache)
    if cache != ""
        @unpack X, b = containers
        data = load(cache)
        parameters.Ts = size(data["X"])[1]
        X[:, 1:parameters.Ts] = data["X"]
        b = data["b"]
    end
end

function save_cache(;containers::SimulationContainers, parameters::SimulationParameters, path, title)
    results_directory = "$path/results"
    simulation_results_directory = "$results_directory/$title"
    !isdir(results_directory) && mkdir(results_directory)
    !isdir(simulation_results_directory) && mkdir(simulation_results_directory)
    save("$(simulation_results_directory)/cache_$(parameters.tmax).jld2" , 
        Dict( 
            "X" => containers.X,
            "b" => containers.b
        )
    )
end