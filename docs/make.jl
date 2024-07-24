using Documenter, Literate, BoreholeNetworksSimulator

pages = [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "API" => "api.md"
]

Literate.markdown("./src/tutorial.jl", "./src")
makedocs(
    pages=pages,
    sitename="BoreholeNetworksSimulator.jl",
    repo=Remotes.GitLab("alblaz", "BoreholeNetworksSimulator")
    )
