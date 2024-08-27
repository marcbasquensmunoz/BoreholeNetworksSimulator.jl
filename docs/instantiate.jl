using Pkg

geothermal_registry = "https://github.com/marcbasquensmunoz/geothermal_registry"

project_directory = joinpath(@__DIR__, "..")
documenter_directory = isempty(ARGS) ? @__DIR__() : joinpath(pwd(), ARGS[1])

cd(project_directory) do
    Pkg.Registry.add(RegistrySpec(url = geothermal_registry))
    Pkg.Registry.add("General")
    Pkg.activate(project_directory)
    Pkg.develop(PackageSpec(path = documenter_directory))
    Pkg.instantiate()
end