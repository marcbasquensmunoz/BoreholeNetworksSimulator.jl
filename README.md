[![documentation (placeholder)](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcbasquensmunoz.github.io/BoreholeNetworksSimulator.jl/stable/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/marcbasquensmunoz/BoreholeNetworksSimulator.jl/actions/workflows/testing.yml/badge.svg?branch=main)](https://github.com/marcbasquensmunoz/BoreholeNetworksSimulator.jl/actions/workflows/testing.yml?query=branch%3Amain)

# BoreholeNetworksSimulator.jl

BoreholeNetworksSimulator.jl is a pure [Julia](https://julialang.org/), performant, and modular framework for simulations of fields of interconnected borehole heat exchangers.
Features:
- Computes fluid temperatures, borehole wall temperatures and heat extracted. 
- Supports many different configurations and settings by being highly modular: boreholes, hydraulic configurations, ground properties, load demand or temperature  constraints, ground boundary conditions can be seamlessly used with minimal effort. 
- Allows for design of operation strategies via an operator callback at each time step.
- Implements the "non-history" time superposition method, reducing the computational complexity in the number of time steps to linear. This allows for simulations with fine time steps.
- Python interoperability.

More information and an extensive list of features can be found in the [documentation](https://marcbasquensmunoz.github.io/BoreholeNetworksSimulator.jl/dev/).


# Quickstart

BoreholeNetworksSimulator.jl is currently not in Julia's General registry (as well as some of its dependencies), however, it is easily available through the local registry [geothermal_registry](https://github.com/marcbasquensmunoz/geothermal_registry). 
In order to install it, start Julia and run the command:

````
using Pkg; pkg"registry add https://github.com/marcbasquensmunoz/geothermal_registry"; Pkg.add("BoreholeNetworksSimulator")
````
