using Documenter, Literate, BoreholeNetworksSimulator

pages = [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "API" => "api.md"
]

dir = @__DIR__
Literate.markdown("$dir/src/tutorial.jl", "$dir/src")
makedocs(
    pages=pages,
    sitename="BoreholeNetworksSimulator.jl",
    #repo=Remotes.GitHub("marcbasquensmunoz", "BoreholeNetworksSimulator.jl")
    repo=Remotes.GitLab("alblaz", "BoreholeNetworksSimulator")
    )
