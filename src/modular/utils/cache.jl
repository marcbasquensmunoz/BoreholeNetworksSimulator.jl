"""
    load_cache!(;containers::SimulationContainers, options::SimulationOptions, cache)

Retrieve the cache stored in `cache` and load it in `containers`. 
Note that the last time step of the previous simulation is loaded in `options.Ts`.
"""
function load_cache!(;containers::SimulationContainers, options::SimulationOptions, cache)
    @unpack Ts = options
    if cache != ""
        @unpack X, b = containers
        data = load(cache)
        Ts = size(data["X"])[1]
        X[:, 1:Ts] = data["X"]
        b = data["b"]
    end
end

"""
    save_cache(;containers::SimulationContainers, path, title)

Save the current simulation results contained in `containers` in a file named `title` in the directory `path`.
"""
function save_cache(;containers::SimulationContainers, path, title)
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