using Documenter, Literate, BoreholeNetworksSimulator

pages = [
    "Introduction" => "index.md"
    #="Tutorial" => [
        "Basic tutorial" => "tutorial.md",
        "Non-history method" => "nonhistory.md",
        "Running from python" => "python.md"
    ],
    "API" => "api.md"=#
]

dir = @__DIR__
#Literate.markdown("$dir/src/tutorial.jl", "$dir/src")
#Literate.markdown("$dir/src/nonhistory.jl", "$dir/src")
makedocs(
    pages = pages,
    sitename = "BoreholeNetworksSimulator.jl"
    #repo = Remotes.GitHub("marcbasquensmunoz", "BoreholeNetworksSimulator.jl")
)
deploydocs(
    repo = Remotes.GitHub("marcbasquensmunoz", "BoreholeNetworksSimulator.jl")
)