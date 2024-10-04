using Pkg

geothermal_registry = "https://github.com/marcbasquensmunoz/geothermal_registry"

project_directory = joinpath(@__DIR__, "..", "..")
plots_directory = joinpath(project_directory, "BNSPlots")

cd(project_directory) do
    Pkg.Registry.add(RegistrySpec(url = geothermal_registry))
    Pkg.Registry.add("General")
    Pkg.activate(; temp=true)
    Pkg.develop(PackageSpec(path = plots_directory))
    Pkg.develop(PackageSpec(path = project_directory))
    Pkg.instantiate()
end