using Pkg

geothermal_registry = "https://github.com/marcbasquensmunoz/geothermal_registry"

documenter_directory = joinpath(@__DIR__, "..")
project_directory = isempty(ARGS) ? @__DIR__() : joinpath(pwd(), ARGS[1])

@show documenter_directory
@show project_directory

cd(project_directory) do
    Pkg.Registry.add(RegistrySpec(url = geothermal_registry))
    Pkg.Registry.add("General")
    Pkg.activate(project_directory)
    Pkg.develop(PackageSpec(path = documenter_directory))
    Pkg.instantiate()
end