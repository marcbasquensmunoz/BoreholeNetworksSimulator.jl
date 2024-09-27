cd(@__DIR__)
import Pkg; Pkg.activate(@__DIR__); Pkg.update(); Pkg.instantiate()
using Documenter, Literate, BoreholeNetworksSimulator, BNSPlots

pages = [
    "Introduction" => "index.md",
    "Tutorial" => [
        "Basic tutorial" => "tutorial.md",
        "Non-history method" => "nonhistory.md",
        "Running from python" => "python.md"
    ],
    "Visualizing the result with BNSPlots" => "BNSPlots.md",
    "Examples" => [
        "Braedstrup borefield" => "Braedstrup.md" 
    ],
    "API" => "api.md"
]

dir = @__DIR__
Literate.markdown("$dir/src/tutorial.jl", "$dir/src")
Literate.markdown("$dir/src/nonhistory.jl", "$dir/src")
makedocs(
    pages = pages,
    sitename = "BoreholeNetworksSimulator.jl"
)
deploydocs(
    repo = "github.com/marcbasquensmunoz/BoreholeNetworksSimulator.jl.git"
)