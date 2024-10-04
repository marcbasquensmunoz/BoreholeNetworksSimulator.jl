using Pkg
Pkg.add(["CSV", "Colors", "Parameters", "WGLMakie"])

project_dir = dirname(pwd())
include("$project_dir/examples/Braedstrup/main.jl")
